from __future__ import annotations

from labw_utils.bioutils.algorithm.sequence import complement, reverse_complement
from labw_utils.bioutils.datastructure.gv import SequenceFuncType, generate_unknown_transcript_id, \
    CanTranscribeInterface
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.record.feature import Feature
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Optional

lh = get_logger(__name__)


class Exon(BaseFeatureProxy, CanTranscribeInterface):
    __slots__ = (
        "_cdna",
    )
    _cdna: Optional[str]

    @property
    def transcript_id(self) -> str:
        return self.attribute_get("transcript_id")


    @property
    def transcribed_length(self):
        return self.naive_length

    def __init__(
            self,
            *,
            data: Feature,
            is_checked: bool,
            shortcut: bool
    ):
        self._cdna = None
        if not shortcut:
            if data.attribute_get("transcript_id") is None:
                data = data.update_attribute(transcript_id=generate_unknown_transcript_id())
        BaseFeatureProxy.__init__(self, data=data, is_checked=is_checked)

    def __repr__(self):
        return f"Exon ({self.start, self.end}) of Transcript {self.transcript_id}"

    def transcribe(self, sequence_func: SequenceFuncType) -> str:
        if self._cdna is None:
            try:
                self._cdna = sequence_func(self.seqname, self.start0b, self.end0b)
                if len(self._cdna) != self.transcribed_length:
                    lh.warning(
                        f"{self.transcript_id}: Different exon length at {self}: " +
                        f"cdna ({len(self._cdna)}) != exon ({self.transcribed_length})"
                    )
                if self.strand is False:
                    self._cdna = reverse_complement(self._cdna)
            except Exception as e:
                lh.warning(f"{self.transcript_id}: Failed to get cDNA sequence at exon ({self.start, self.end}) {e}")
                self._cdna = ""
        return self._cdna
