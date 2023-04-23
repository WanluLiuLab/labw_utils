from __future__ import annotations

import bisect
import math
import operator

from labw_utils.bioutils.datastructure.gv import DEFAULT_SORT_EXON_EXON_STRAND_POLICY, GVPError, CanTranscribeInterface, \
    SortedContainerInterface
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy, update_gene_id, update_transcript_id
from labw_utils.bioutils.record.feature import FeatureInterface, FeatureType
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Optional, Iterable, Tuple, Union, Callable, SequenceProxy

lh = get_logger(__name__)


class ExonInATranscriptOnDifferentChromosomeError(GVPError):
    pass


class ExonInATranscriptOnDifferentStrandError(GVPError):
    pass


class Transcript(
    BaseFeatureProxy,
    CanTranscribeInterface,
    SortedContainerInterface
):
    """
    Transcript is a list of exons, always sorted.
    """

    __slots__ = [
        "_exons",
        "_cdna",
        "_is_inferred",
        "_transcript_id",
        "_gene_id"
    ]
    _exons: List[Exon]
    _cdna: Optional[str]
    _is_inferred: bool
    _exon_boundaries: Optional[List[Tuple[int, int]]]
    _splice_sites: Optional[List[Tuple[int, int]]]
    _transcript_id: str
    _gene_id: str

    @property
    def transcript_id(self) -> str:
        return self._transcript_id

    @property
    def gene_id(self) -> str:
        return self._gene_id

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
    def exon_boundaries(self) -> SequenceProxy[Tuple[int, int]]:
        if self._exon_boundaries is None:
            self._exon_boundaries = list((exon.start, exon.end) for exon in self._exons)
        return SequenceProxy(self._exon_boundaries)

    @property
    def splice_sites(self) -> SequenceProxy[Tuple[int, int]]:
        if self._splice_sites is None:
            self._splice_sites = list(
                (self._exons[i].end, self._exons[i + 1].start)
                for i in range(self.number_of_exons - 1)
            )
        return SequenceProxy(self._splice_sites)

    @property
    def exons(self) -> SequenceProxy[Exon]:
        """Get Exon Iterator"""
        return SequenceProxy(self._exons)

    def __init__(
            self,
            *,
            data: FeatureInterface,
            is_checked: bool,
            keep_sorted: bool,
            shortcut: bool,
            exons: Iterable[Exon],
            is_inferred: bool
    ):
        self._is_sorted = keep_sorted
        self._cdna = None
        self._exon_boundaries = None
        self._splice_sites = None
        self._is_inferred = is_inferred
        if not shortcut:
            self._transcript_id, data = update_transcript_id(data)
            self._gene_id, data = update_gene_id(data)
            if self._is_inferred and data.parsed_feature is not FeatureType.TRANSCRIPT:
                data = data.update(feature="transcript")
            self._exons = list(exons)
        else:
            self._transcript_id = data.attribute_get("transcript_id")  # type: ignore
            self._gene_id = data.attribute_get("gene_id")  # type: ignore
            self._exons = exons  # type: ignore
        BaseFeatureProxy.__init__(self, data=data, is_checked=is_checked)

    def __repr__(self):
        return f"Transcript {self.transcript_id} of gene {self.gene_id}"

    def exon_level_equiv(self, other: Transcript) -> bool:
        if self.number_of_exons != other.number_of_exons:
            return False
        return all(map(lambda exon_pair: operator.eq(*exon_pair), zip(self.exons, other.exons)))

    def get_exon(self, exon_id: int) -> Exon:
        return self._exons[exon_id]

    def get_intron_length(self, intron_index: int) -> Union[int, float]:
        if intron_index == -1 or intron_index == self.number_of_exons:
            return math.inf
        return - operator.sub(*self.exon_boundaries[intron_index])

    def update_exon_number(
            self,
            exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY
    ) -> Transcript:
        """
        Renew exon numbers.

        :param exon_number_policy: The UCSC style. TODO
        """

        def update_exon_number(_exon: Exon, target_exon_number: int) -> Exon:
            return Exon(
                data=_exon.get_data().update_attribute(exon_number=target_exon_number),
                is_checked=_exon.is_checked,
                shortcut=True
            )

        new_exons = []
        if exon_number_policy == "stranded":
            if self.strand is True:
                for i in range(self.number_of_exons):
                    new_exons.append(update_exon_number(self._exons[i], i + 1))
            elif self.strand is False:
                for i in range(self.number_of_exons):
                    new_exons.append(
                        update_exon_number(self._exons[self.number_of_exons - i - 1], i + 1)
                    )
                new_exons.reverse()
            else:
                raise ValueError("Unstranded exon detected!")
        elif exon_number_policy == "unstranded":
            for i in range(self.number_of_exons):
                new_exons.append(update_exon_number(self._exons[i], i + 1))
        return Transcript(
            data=self._data,
            is_checked=self._is_checked,
            keep_sorted=self._is_sorted,
            exons=new_exons,
            shortcut=True,
            is_inferred=False
        )

    def rescale_from_exon_boundaries(self, force: bool = False) -> Transcript:
        """
        This method is only used if the transcript is inferred from exon.

        For example, in a GTF that contains only exons,
        we need to infer transcript and gene from existing exons.

        :param force: Force to rescale transcript from exon boundaries.
        :returns: New instance whose boundary is rescaled.
        """
        if self._is_inferred or force:
            new_data = self._data.update(start=self._exons[0].start, end=self._exons[-1].end)
            return Transcript(
                data=new_data,
                is_checked=self._is_checked,
                keep_sorted=self._is_sorted,
                exons=self._exons,
                shortcut=True,
                is_inferred=False
            )
        else:
            return self

    def add_exon(self, exon: Exon) -> Transcript:
        new_exons = list(self._exons)
        if not self._is_checked:
            if exon.seqname != self.seqname:
                raise ExonInATranscriptOnDifferentChromosomeError
            if exon.strand != self.strand and exon.strand is not None:
                raise ExonInATranscriptOnDifferentStrandError
        if self._is_sorted:
            new_pos = bisect.bisect_left(new_exons, exon)
            new_exons.insert(new_pos, exon)
        else:
            new_exons.append(exon)
        return Transcript(
            data=self._data,
            is_checked=self._is_checked,
            keep_sorted=self._is_sorted,
            exons=new_exons,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def del_exon(self, exon_index: int) -> Transcript:
        """
        Delete an exon with corresponding index.

        .. warning::
            Exon index is NOT exon number!

            * If the transcript is positively stranded, checked and sorted with exon number rearranged,
              Exon index is Exon number - 1.
            * If the transcript is negatively stranded, checked and sorted with exon number rearranged,
              Exon index is Exon number - 1 in reversed order.
            * Unstranded similar to positively stranded.
        """
        new_exons = list(self._exons)
        _ = new_exons.pop(exon_index)
        return Transcript(
            data=self._data,
            is_checked=self._is_checked,
            keep_sorted=self._is_sorted,
            exons=new_exons,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def transcribe(self, sequence_func: Callable[[str, int, int], str]) -> str:
        if self._cdna is None:
            if self.strand is False:
                self._cdna = "".join(exon.transcribe(sequence_func) for exon in self._exons[::-1])
            else:
                self._cdna = "".join(exon.transcribe(sequence_func) for exon in self._exons)
            if len(self._cdna) != self.transcribed_length:
                lh.warning(
                    f"Transcript {self.transcript_id} " +
                    f"cdna_len({len(self._cdna)}) != transcribed_len ({self.transcribed_length})."
                )
        return self._cdna
