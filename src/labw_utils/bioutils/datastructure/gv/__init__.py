from __future__ import annotations

import uuid
from typing import Callable, TypeVar

VALID_SORT_EXON_EXON_STRAND_POLICY = ("unstranded", "stranded", "none")
DEFAULT_SORT_EXON_EXON_STRAND_POLICY = "unstranded"
SequenceFuncType = Callable[[str, int, int], str]


class _NotSet:
    pass


_notset = _NotSet()
_T = TypeVar("_T")


def generate_unknown_transcript_id() -> str:
    """Generate a new unknown transcript ID"""
    return 'unknown_transcript_id' + str(uuid.uuid4())


def generate_unknown_gene_id() -> str:
    """Generate a new unknown gene ID"""
    return 'unknown_gene_id' + str(uuid.uuid4())


class GVPError(ValueError):
    pass


class CanTranscribe:
    def transcribe(self, sequence_func: SequenceFuncType) -> str:
        raise NotImplementedError
