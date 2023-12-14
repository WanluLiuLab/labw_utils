"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

import bisect
import itertools

from labw_utils.bioutils.algorithm.utils import merge_intervals
from labw_utils.bioutils.datastructure.gv import GVPError
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy, update_gene_id
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container_interface import (
    TranscriptContainerInterface,
    DuplicatedTranscriptIDError,
)
from labw_utils.bioutils.record.feature import FeatureType, FeatureInterface
from labw_utils.typing_importer import List, Optional, Iterable, Sequence

from labw_utils.typing_importer import SequenceProxy


class TranscriptInAGeneOnDifferentChromosomeError(GVPError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    def __init__(self, transcript: Transcript, gene_seqname: str):
        super().__init__(
            f"{transcript}: " f"gene.seqname={gene_seqname}, " f"while transcript.seqname={transcript.seqname}"
        )


class TranscriptInAGeneOnDifferentStrandError(GVPError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    def __init__(self, transcript: Transcript, gene_strand: Optional[bool]):
        super().__init__(
            f"{transcript}: " f"gene.strand={gene_strand}, " f"while transcript.strand={transcript.strand}"
        )


class Gene(BaseFeatureProxy, TranscriptContainerInterface):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    __slots__ = ("_transcripts", "_transcript_ids", "_gene_id")
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

    @property
    def transcribed_length(self) -> int:
        return sum(transcript.transcribed_length for transcript in self._transcripts)

    @property
    def mappable_length(self) -> int:
        all_exons: List[Exon] = list(itertools.chain(*list(transcript.exons for transcript in self._transcripts)))
        all_intervals = list((exon.start, exon.end) for exon in all_exons)
        merged_intervals = merge_intervals(all_intervals)
        return sum(interval[1] - interval[0] + 1 for interval in merged_intervals)

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self._transcripts[self._transcript_ids.index(transcript_id)]

    def __init__(
        self,
        *,
        data: FeatureInterface,
        is_checked: bool,
        shortcut: bool,
        transcripts: Optional[Iterable[Transcript]],
        transcript_ids: Optional[Iterable[str]],
        is_inferred: bool,
    ):
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
        BaseFeatureProxy.__init__(
            self,
            data=data,
            is_checked=self._is_inferred,
            shortcut=shortcut,
            transcripts=self._transcripts,
            transcript_ids=self._transcript_ids,
            is_inferred=self._is_inferred,
        )

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
        new_pos = bisect.bisect_left(self._transcripts, transcript)
        new_transcripts.insert(new_pos, transcript)
        new_transcript_ids.insert(new_pos, transcript.transcript_id)
        return Gene(
            data=self._data,
            is_checked=self._is_checked,
            is_inferred=self._is_inferred,
            transcripts=new_transcripts,
            transcript_ids=new_transcript_ids,
            shortcut=True,
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
            transcripts=new_transcripts,
            transcript_ids=new_transcript_ids,
            is_inferred=self._is_inferred,
            shortcut=True,
        )

    def replace_transcript(self, new_transcript: Transcript) -> Gene:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)

    def __repr__(self):
        return f"Gene {self.gene_id}"

    def collapse_transcript(
        self, skip_isoform_on_different_contig_or_strand_behaviour: bool = False
    ) -> Optional[Transcript]:
        transcripts = []
        contig_of_first_transcript = self.seqname
        strand_of_first_transcript = self.strand
        for transcript in self._transcripts:
            if len(transcript.exons) == 0:
                continue
            if contig_of_first_transcript != transcript.seqname:
                if skip_isoform_on_different_contig_or_strand_behaviour:
                    continue
                else:
                    raise ValueError
            if strand_of_first_transcript != transcript.strand:
                if skip_isoform_on_different_contig_or_strand_behaviour:
                    continue
                else:
                    raise ValueError
            transcripts.append(transcript)
        all_exon_bondary = merge_intervals(
            list(itertools.chain(*(list(transcript.exon_boundaries) for transcript in transcripts)))
        )
        exon_list = [
            Exon(
                data=(
                    self._data.update(
                        feature="exon",
                        strand=strand_of_first_transcript,
                        seqname=contig_of_first_transcript,
                        start=exon_bondary[0] + 1,
                        end=exon_bondary[1],
                    ).update_attribute(transcript_id=self.gene_id)
                ),
                is_checked=self.is_checked,
                shortcut=True,
            )
            for exon_bondary in all_exon_bondary
        ]
        if len(exon_list) == 0:
            exon_list = [
                Exon(
                    data=(
                        self._data.update(
                            feature="exon",
                            strand=strand_of_first_transcript,
                            seqname=contig_of_first_transcript,
                        ).update_attribute(transcript_id=self.gene_id)
                    ),
                    is_checked=self.is_checked,
                    shortcut=True,
                )
            ]
        return Transcript(
            data=(
                self._data.update(
                    feature="transcript",
                    strand=strand_of_first_transcript,
                    seqname=contig_of_first_transcript,
                ).update_attribute(transcript_id=self.gene_id)
            ),
            is_checked=self.is_checked,
            is_inferred=True,
            shortcut=True,
            exons=sorted(exon_list, key=lambda e: e.start0b),
        ).rescale_from_exon_boundaries()

    def gc(self):
        for transcript in self._transcripts:
            transcript.gc()


class DumbGene(Gene):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    def __init__(
        self,
        *,
        data: FeatureInterface,
        is_checked: bool,
        shortcut: bool,
        transcripts: Optional[Iterable[Transcript]],
        transcript_ids: Optional[Iterable[str]],
        is_inferred: bool,
    ):
        _ = is_checked
        del is_checked

        Gene.__init__(
            self,
            data=data,
            is_checked=True,
            shortcut=shortcut,
            transcripts=transcripts,
            transcript_ids=transcript_ids,
            is_inferred=is_inferred,
        )

    def __repr__(self):
        return f"DumbGene {self.gene_id}"
