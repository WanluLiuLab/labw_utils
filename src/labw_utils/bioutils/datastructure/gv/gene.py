from __future__ import annotations

import bisect

from labw_utils.bioutils.datastructure.gv import GVPError, SortedContainerInterface
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy, update_gene_id
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container_interface import TranscriptContainerInterface, \
    DuplicatedTranscriptIDError
from labw_utils.bioutils.record.feature import FeatureType, FeatureInterface
from labw_utils.typing_importer import List, Optional, Iterable, Sequence, SequenceProxy


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


class Gene(BaseFeatureProxy, TranscriptContainerInterface, SortedContainerInterface):
    __slots__ = (
        "_transcripts",
        "_transcript_ids",
        "_gene_id"
    )
    _gene_id: str
    _transcripts: List[Transcript]
    _transcript_ids: List[str]

    @property
    def number_of_transcripts(self) -> int:
        return len(self._transcripts)

    @property
    def transcript_values(self) -> Sequence[Transcript]:
        return SequenceProxy(self._transcripts)

    @property
    def transcript_ids(self) -> Sequence[str]:
        return SequenceProxy(self._transcript_ids)

    @property
    def gene_id(self) -> str:
        return self._gene_id

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self._transcripts[self._transcript_ids.index(transcript_id)]

    def __init__(
            self,
            *,
            data: FeatureInterface,
            is_checked: bool,
            keep_sorted: bool,
            shortcut: bool,
            transcripts: Optional[Iterable[Transcript]],
            transcript_ids: Optional[Iterable[str]],
            is_inferred: bool
    ):
        self._is_sorted = keep_sorted
        self._is_inferred = is_inferred
        if transcripts is None:
            transcripts = []
        if transcript_ids is None:
            transcript_ids = []
        if not shortcut:
            self._gene_id, data = update_gene_id(data)
            if self._is_inferred and data.parsed_feature is not FeatureType.GENE:
                data = data.update(feature="gene")
            self._transcripts = list(transcripts)
            self._transcript_ids = list(transcript_ids)
        else:
            self._gene_id = data.attribute_get("gene_id")  # type: ignore
            self._transcripts = transcripts  # type: ignore
            self._transcript_ids = transcript_ids  # type: ignore
        BaseFeatureProxy.__init__(self, data=data, is_checked=is_checked)

    def add_transcript(self, transcript: Transcript) -> Gene:
        if not self._is_checked:
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


class DumbGene(Gene):

    def __init__(
            self,
            *,
            data: FeatureInterface,
            is_checked: bool,
            keep_sorted: bool,
            shortcut: bool,
            transcripts: Optional[Iterable[Transcript]],
            transcript_ids: Optional[Iterable[str]],
            is_inferred: bool
    ):
        _ = is_checked, keep_sorted
        del is_checked, keep_sorted

        Gene.__init__(
            self,
            data=data,
            is_checked=True,
            keep_sorted=False,
            shortcut=shortcut,
            transcripts=transcripts,
            transcript_ids=transcript_ids,
            is_inferred=is_inferred
        )

    def __repr__(self):
        return f"DumbGene {self.gene_id}"
