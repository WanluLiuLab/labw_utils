from __future__ import annotations

from typing import Final, List, Optional, Iterable

from labw_utils.bioutils.datastructure.gv import generate_unknown_gene_id, GVPError
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.gv.transcript_container import TranscriptContainer, TranscriptContainerInterface
from labw_utils.bioutils.record.feature import Feature, FeatureType


class TranscriptInAGeneOnDifferentChromosomeError(GVPError):
    pass


class DuplicatedTranscriptError(GVPError):
    pass


class TranscriptInAGeneOnDifferentStrandError(GVPError):
    pass


class Gene(BaseFeatureProxy, TranscriptContainerInterface):
    _transcript_container: TranscriptContainer
    preserved_attributes: Final[List[str]] = ("gene_id",)

    @property
    def number_of_transcripts(self) -> int:
        return self._transcript_container.number_of_transcripts

    @property
    def transcript_values(self) -> Iterable[Transcript]:
        return self._transcript_container.transcript_values

    @property
    def transcript_ids(self) -> Iterable[str]:
        return self._transcript_container.transcript_ids

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self._transcript_container.get_transcript(transcript_id)

    @property
    def gene_id(self) -> str:
        return self.attribute_get("gene_id")

    def __init__(
            self,
            data: Feature,
            shortcut: bool = False,
            transcript_container: Optional[TranscriptContainer] = None,
            is_inferred: bool = False,
            **kwargs
    ):
        if transcript_container is None:
            transcript_container = TranscriptContainer([])
        self._is_inferred = is_inferred
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
        self._transcript_container = transcript_container

    def add_transcript(self, transcript: Transcript) -> Gene:
        if transcript.seqname != self.seqname:
            raise TranscriptInAGeneOnDifferentChromosomeError(
                f"gene.seqname={self.seqname}, while transcript.seqname={transcript.seqname}"
            )
        if transcript.strand != self.strand and transcript.strand is not None:
            raise TranscriptInAGeneOnDifferentStrandError
        new_transcript_container = self._transcript_container.add_transcript(transcript)
        return Gene(
            data=self._data,
            transcript_container=new_transcript_container,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def del_transcript(self, transcript_id: str) -> Gene:
        new_transcript_container = self._transcript_container.del_transcript(transcript_id)
        return Gene(
            data=self._data,
            transcript_container=new_transcript_container,
            is_inferred=self._is_inferred,
            shortcut=True
        )

    def replace_transcript(self, new_transcript: Transcript) -> Gene:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)

    def __repr__(self):
        return f"Gene {self.gene_id}"
