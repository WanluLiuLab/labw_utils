from __future__ import annotations

from typing import Optional, Final, List

from labw_utils.bioutils.datastructure.gv import SequenceFuncType, generate_unknown_transcript_id, \
    generate_unknown_gene_id, CanTranscribe
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.record.feature import Feature
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


class Exon(BaseFeatureProxy, CanTranscribe):
    __slots__ = (
        "_cdna",
    )
    _cdna: Optional[str]
    _preserved_attrs: Final[List[str]] = ("gene_id", "transcript_id", "exon_number")

    @property
    def transcript_id(self) -> str:
        return self.attribute_get("transcript_id")

    @property
    def gene_id(self) -> str:
        return self.attribute_get("gene_id")

    @property
    def exon_number(self) -> int:
        return self.attribute_get("exon_number")

    @property
    def transcribed_length(self):
        return self.naive_length

    def __init__(self, data: Feature, shortcut: bool = False, **kwargs):
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
            if data.attribute_get("exon_number") is None:
                should_update_attributes = True
                data_attributes_update_kwargs["exon_number"] = 0
            if should_update_attributes:
                data = data.update(attribute=data_attributes_update_kwargs)
        super().__init__(data, **kwargs)
        self._cdna = None

    def __repr__(self):
        return f"Exon {self.exon_number} of {self.transcript_id}"

    def transcribe(self, sequence_func: SequenceFuncType) -> str:
        if self._cdna is None:
            try:
                self._cdna = sequence_func(self.seqname, self.start0b, self.end0b)
                if len(self._cdna) != self.transcribed_length:
                    lh.warn(
                        f"{self.transcript_id}: Different exon length at {self}: " +
                        f"cdna ({len(self._cdna)}) != exon ({self.transcribed_length})"
                    )
            except Exception as e:  # TODO
                lh.warn(f"{self.transcript_id}: Failed to get cDNA sequence at exon {self.exon_number} {e}")
                self._cdna = ""
        return self._cdna
