"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

from labw_utils.bioutils.algorithm.sequence import reverse_complement
from labw_utils.bioutils.datastructure.gv import (
    SequenceFuncType,
    CanTranscribeInterface,
    dumb_legalize_region_func,
    LegalizeRegionFuncType,
)
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy, update_transcript_id
from labw_utils.bioutils.record.feature import FeatureInterface
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Optional

_lh = get_logger(__name__)


class Exon(BaseFeatureProxy, CanTranscribeInterface):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    __slots__ = ("_cdna", "_transcript_id")
    _cdna: Optional[str]
    _transcript_id: str

    @property
    def transcript_id(self) -> str:
        return self._transcript_id

    @property
    def transcribed_length(self):
        return self.naive_length

    def __init__(self, *, data: FeatureInterface, is_checked: bool, shortcut: bool):
        self._cdna = None
        if not shortcut:
            self._transcript_id, data = update_transcript_id(data)
        else:
            self._transcript_id = data.attribute_get("transcript_id")  # type: ignore
        BaseFeatureProxy.__init__(self, data=data, is_checked=is_checked)

    def __repr__(self):
        return f"Exon ({self.start, self.end}) of Transcript {self.transcript_id}"

    def transcribe(
        self, sequence_func: SequenceFuncType, legalize_region_func: LegalizeRegionFuncType = dumb_legalize_region_func
    ) -> str:
        if self._cdna is None:
            try:
                self._cdna = sequence_func(*legalize_region_func(self.seqname, self.start0b, self.end0b))
                if len(self._cdna) != self.transcribed_length:
                    _lh.warning(
                        f"{self.transcript_id}: Different exon length at {self}: "
                        + f"cdna ({len(self._cdna)}) != exon ({self.transcribed_length})"
                    )
                if self.strand is False:
                    self._cdna = reverse_complement(self._cdna)
            except Exception as e:
                _lh.warning(f"{self.transcript_id}: Failed to get cDNA sequence at exon ({self.start, self.end}) {e}")
                self._cdna = ""
        return self._cdna

    def gc(self):
        self._cdna = None
