from __future__ import annotations

import bisect
from typing import Final, List, Optional, Iterable

from labw_utils.bioutils.datastructure.gv import generate_unknown_gene_id, GVPError, ContainerInterface
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container import TranscriptContainerInterface
from labw_utils.bioutils.record.feature import Feature, FeatureType


class TranscriptInAGeneOnDifferentChromosomeError(GVPError):
    pass


class DuplicatedTranscriptError(GVPError):
    pass


class TranscriptInAGeneOnDifferentStrandError(GVPError):
    pass


class DuplicatedTranscriptIDError(GVPError):
    pass


class Gene(BaseFeatureProxy, TranscriptContainerInterface, ContainerInterface):
    _gene_id: str
    preserved_attributes: Final[List[str]] = ("gene_id",)
    __slots__ = (
        "_transcripts",
        "_transcript_ids"
    )
    _transcripts: List[Transcript]
    _transcript_ids: List[str]

    @property
    def number_of_transcripts(self) -> int:
        return len(self._transcripts)

    @property
    def transcript_values(self) -> Iterable[Transcript]:
        return iter(self._transcripts)

    @property
    def transcript_ids(self) -> Iterable[str]:
        return iter(self._transcript_ids)

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self._transcripts[self._transcript_ids.index(transcript_id)]

    @property
    def gene_id(self) -> str:
        return self._gene_id

    def __init__(
            self,
            data: Feature,
            *,
            keep_sorted: bool = False,
            shortcut: bool = False,
            transcripts: Optional[Iterable[Transcript]] = None,
            transcript_ids: Optional[Iterable[str]] = None,
            is_inferred: bool = False,
            **kwargs
    ):
        self._is_sorted = keep_sorted
        self._is_inferred = is_inferred
        if transcripts is None:
            self._transcripts = []
        else:
            self._transcripts = list(transcripts)
        if transcript_ids is None:
            self._transcript_ids = list(transcript.transcript_id for transcript in self._transcripts)
        else:
            self._transcript_ids = list(transcript_ids)
        if not shortcut:
            should_update_attributes = False
            data_attributes_update_kwargs = {k: v for k, v in zip(
                data.attribute_keys, data.attribute_values
            )}
            if data.attribute_get("gene_id") is None:
                should_update_attributes = True
                data_attributes_update_kwargs["gene_id"] = generate_unknown_gene_id()
            if self._is_inferred and data.parsed_feature is not FeatureType.Gene:
                data = data.update(feature="gene")
            if should_update_attributes:
                data = data.update(attribute=data_attributes_update_kwargs)
        BaseFeatureProxy.__init__(self, data, **kwargs)
        self._gene_id = self.attribute_get("gene_id")

    def add_transcript(self, transcript: Transcript) -> Gene:
        if transcript.seqname != self.seqname:
            raise TranscriptInAGeneOnDifferentChromosomeError(
                f"gene.seqname={self.seqname}, while transcript.seqname={transcript.seqname}"
            )
        if transcript.strand != self.strand and transcript.strand is not None:
            raise TranscriptInAGeneOnDifferentStrandError
        new_transcripts = list(self._transcripts)
        new_transcript_ids = list(self._transcript_ids)
        if transcript.transcript_id in self._transcript_ids:
            raise DuplicatedTranscriptIDError(
                f"Transcript ID {transcript.transcript_id} duplicated"
            )
        if self._is_sorted:
            new_pos = bisect.bisect_left(self._transcripts, transcript)
            new_transcripts.insert(new_pos, transcript)
            new_transcript_ids.insert(new_pos, transcript.transcript_id)
        else:
            new_transcripts.append(transcript)
            new_transcript_ids.append(transcript.transcript_id)

        return Gene(
            data=self._data,
            keep_sorted=self._is_sorted,
            transcripts=new_transcripts,
            transcript_ids=self.transcript_ids,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def del_transcript(self, transcript_id: str) -> Gene:
        new_transcripts = list(self._transcripts)
        new_transcript_ids = list(self._transcript_ids)
        pop_pos = self._transcript_ids.index(transcript_id)
        _ = new_transcripts.pop(pop_pos)
        _ = new_transcript_ids.pop(pop_pos)
        return Gene(
            data=self._data,
            keep_sorted=self._is_sorted,
            transcripts=new_transcripts,
            transcript_ids=self.transcript_ids,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def replace_transcript(self, new_transcript: Transcript) -> Gene:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)

    def __repr__(self):
        return f"Gene {self.gene_id}"
