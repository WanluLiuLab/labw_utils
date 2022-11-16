from __future__ import annotations

import bisect
import math
from typing import Final, List, Optional, Iterable, Tuple, Union, Callable

from labw_utils.bioutils.algorithm.sequence import reverse_complement
from labw_utils.bioutils.datastructure.gv import DEFAULT_SORT_EXON_EXON_STRAND_POLICY, generate_unknown_transcript_id, \
    generate_unknown_gene_id, GVPError, CanTranscribeInterface
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.record.feature import Feature, FeatureType
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


class ExonInATranscriptOnDifferentChromosomeError(GVPError):
    pass


class DuplicatedExonError(GVPError):
    pass


class ExonInATranscriptOnDifferentStrandError(GVPError):
    pass


class Transcript(BaseFeatureProxy, CanTranscribeInterface):
    """
    Transcript is a list of exons, always sorted.
    """
    preserved_attributes: Final[List[str]] = ("gene_id", "transcript_id")
    __slots__ = (
        "_exons",
        "_cdna",
        "_is_inferred"
    )
    _exons: List[Exon]
    _cdna: Optional[str]
    _is_inferred: Optional[bool]

    @property
    def transcript_id(self) -> str:
        return self.attribute_get("transcript_id")

    @property
    def gene_id(self) -> str:
        return self.attribute_get("gene_id")

    @property
    def number_of_exons(self) -> int:
        return len(self._exons)

    @property
    def span_length(self) -> int:
        """
        The spanning length of all exons
        """
        if self.number_of_exons == 0:
            return -1
        else:
            return self._exons[-1].end0b - self._exons[0].start0b

    @property
    def transcribed_length(self) -> int:
        """
        Length after transcribed to cDNA
        """
        return sum(exon.transcribed_length for exon in self._exons)

    @property
    def exon_boundaries(self) -> Iterable[Tuple[int, int]]:
        for exon in self._exons:
            yield exon.start, exon.end

    @property
    def splice_sites(self) -> Iterable[Tuple[int, int]]:
        for i in range(len(self._exons) - 1):
            yield self._exons[i].end, self._exons[i + 1].start

    @property
    def exons(self) -> Iterable[Exon]:
        """Get Exon Iterator"""
        return iter(self._exons)

    def __init__(
            self,
            data: Feature,
            shortcut: bool = False,
            exons: Optional[List[Exon]] = None,
            is_inferred: bool = False,
            **kwargs
    ):
        if exons is None:
            exons = []
        self._is_inferred = is_inferred
        if not shortcut:
            should_update_attributes = False
            data_attributes_update_kwargs = {k: v for k, v in zip(
                data.attribute_keys, data.attribute_values
            )}
            if data.attribute_get("transcript_id") is None:
                should_update_attributes = True
                data_attributes_update_kwargs["transcript_id"] = generate_unknown_transcript_id()
            if data.attribute_get("gene_id") is None:
                should_update_attributes = True
                data_attributes_update_kwargs["gene_id"] = generate_unknown_gene_id()
            if self._is_inferred and data.parsed_feature is not FeatureType.Transcript:
                data = data.update(feature="transcript")
            if should_update_attributes:
                data = data.update(attribute=data_attributes_update_kwargs)
        super().__init__(data, **kwargs)
        self._cdna = None
        self._exons = list(exons)

    def __repr__(self):
        return f"Transcript {self.transcript_id} of {self.gene_id}"

    def exon_level_equiv(self, other: Transcript) -> bool:
        for exon_s, exon_o in zip(self._exons, other._exons):
            if not exon_s == exon_o:
                return False
        return True

    def get_exon(self, exon_id: int) -> Exon:
        return self._exons[exon_id]

    def get_intron_length(self, intron_index: int) -> Union[int, float]:
        if intron_index == -1 or intron_index == self.number_of_exons:
            return math.inf
        _len = self._exons[intron_index + 1].start - self._exons[intron_index].end + 1
        return _len

    def update_exon_number(
            self,
            exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY
    ) -> Transcript:
        if exon_number_policy == "stranded":
            if self.strand is True:
                for i in range(len(self._exons)):
                    self._exons[i] = self._exons[i].update(exon_number=i + 1)
            elif self.strand is False:
                for i in range(len(self._exons)):
                    self._exons[len(self._exons) - i - 1] = \
                        self._exons[len(self._exons) - i - 1].update(exon_number=i + 1)
            else:
                raise ValueError("Unstranded exon detected!")
        elif exon_number_policy == "unstranded":
            for i in range(len(self._exons)):
                self._exons[i] = self._exons[i].update(exon_number=i + 1)
        return self

    def rescale_from_exon_boundaries(self) -> Transcript:
        if self._is_inferred:
            new_data = self._data.update(start=self._exons[0].start, end=self._exons[-1].end)
            return Transcript(
                data=new_data,
                exons=self._exons,
                shortcut=True,
                is_inferred=False
            )
        else:
            return self

    @classmethod
    def infer_from_exon(cls, exon: Exon) -> Transcript:
        return BaseFeatureProxy.cast_to(exon, class_type=Transcript, is_inferred=True)

    def add_exon(self, exon: Exon) -> Transcript:
        new_exons = list(self._exons)
        new_pos = bisect.bisect_left(new_exons, exon)
        if new_pos < len(new_exons) and exon == new_exons[new_pos]:
            raise DuplicatedExonError
        if exon.seqname != self.seqname:
            raise ExonInATranscriptOnDifferentChromosomeError
        if exon.strand != self.strand and exon.strand is not None:
            raise ExonInATranscriptOnDifferentStrandError
        new_exons.insert(new_pos, exon)
        return Transcript(
            data=self._data,
            exons=new_exons,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def del_exon(self, exon_number: int) -> Transcript:
        new_exons = list(self._exons)
        _ = new_exons.pop(exon_number)
        return Transcript(
            data=self._data,
            exons=new_exons,
            shortcut=True
        )

    def transcribe(self, sequence_func: Callable[[str, int, int], str]) -> str:
        if self._cdna is None:
            if self.strand is False:
                self._cdna = "".join(reverse_complement(exon.transcribe(sequence_func)) for exon in self._exons[::-1])

            else:
                self._cdna = "".join(exon.transcribe(sequence_func) for exon in self._exons)
            if len(self._cdna) != self.transcribed_length:
                lh.warn(
                    f"Transcript {self.transcript_id} " +
                    f"cdna_len({len(self._cdna)}) != transcribed_len ({self.transcribed_length})."
                )
        return self._cdna
