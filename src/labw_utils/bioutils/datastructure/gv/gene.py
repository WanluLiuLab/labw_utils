from __future__ import annotations

import bisect
from typing import Final, List, Optional, Iterable

from labw_utils.bioutils.datastructure.gv import generate_unknown_gene_id, GVPError, SortedContainerInterface
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container_interface import TranscriptContainerInterface
from labw_utils.bioutils.record.feature import Feature, FeatureType


class TranscriptInAGeneOnDifferentChromosomeError(GVPError):
    def __init__(self, transcript: Transcript, gene_seqname: str):
        super().__init__(
            f"{transcript}: "
            f"gene.seqname={gene_seqname}, "
            f"while transcript.seqname={transcript.seqname}"
        )


class TranscriptInAGeneOnDifferentStrandError(GVPError):
    def __init__(self, transcript: Transcript, gene_strand: Optional[bool]):
        super().__init__(
            f"{transcript}: "
            f"gene.strand={gene_strand}, "
            f"while transcript.strand={transcript.strand}"
        )


class DuplicatedTranscriptIDError(GVPError):
    def __init__(self, transcript_id: str):
        super().__init__(
            f"Transcript ID {transcript_id} duplicated"
        )


class Gene(BaseFeatureProxy, TranscriptContainerInterface, SortedContainerInterface):
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

    @property
    def gene_id(self) -> str:
        return self._gene_id

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self._transcripts[self._transcript_ids.index(transcript_id)]

    def __init__(
            self,
            *,
            data: Feature,
            is_checked: bool,
            keep_sorted: bool,
            shortcut: bool,
            transcripts: Optional[Iterable[Transcript]],
            transcript_ids: Optional[Iterable[str]],
            is_inferred: bool
    ):
        self._is_sorted = keep_sorted
        self._is_inferred = is_inferred
        self._transcripts = list(transcripts)
        self._transcript_ids = list(transcript_ids)
        if not shortcut:
            should_update_attributes = False
            data_attributes_update_kwargs = {k: v for k, v in zip(
                data.attribute_keys, data.attribute_values
            )}
            if data.attribute_get("gene_id") is None:
                should_update_attributes = True
                data_attributes_update_kwargs["gene_id"] = generate_unknown_gene_id()
            if self._is_inferred and data.parsed_feature is not FeatureType.GENE:
                data = data.update(feature="gene")
            if should_update_attributes:
                data = data.update(attribute=data_attributes_update_kwargs)
        BaseFeatureProxy.__init__(self, data=data, is_checked=is_checked)
        self._gene_id = self.attribute_get("gene_id")

    def add_transcript(self, transcript: Transcript) -> Gene:
        if self._is_checked:
            if transcript.seqname != self.seqname:
                raise TranscriptInAGeneOnDifferentChromosomeError(transcript, self.seqname)
            if transcript.strand != self.strand and transcript.strand is not None:
                raise TranscriptInAGeneOnDifferentStrandError(transcript, self.strand)
            if transcript.transcript_id in self._transcript_ids:
                raise DuplicatedTranscriptIDError(transcript.transcript_id)
        new_transcripts = list(self._transcripts)
        new_transcript_ids = list(self._transcript_ids)
        if self._is_sorted:
            new_pos = bisect.bisect_left(self._transcripts, transcript)
            new_transcripts.insert(new_pos, transcript)
            new_transcript_ids.insert(new_pos, transcript.transcript_id)
        else:
            new_transcripts.append(transcript)
            new_transcript_ids.append(transcript.transcript_id)
        return Gene(
            data=self._data,
            is_checked=self._is_checked,
            keep_sorted=self._is_sorted,
            is_inferred=self._is_inferred,
            transcripts=new_transcripts,
            transcript_ids=new_transcript_ids,
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
            is_checked=self._is_checked,
            keep_sorted=self._is_sorted,
            transcripts=new_transcripts,
            transcript_ids=new_transcript_ids,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def replace_transcript(self, new_transcript: Transcript) -> Gene:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)

    def __repr__(self):
        return f"Gene {self.gene_id}"
