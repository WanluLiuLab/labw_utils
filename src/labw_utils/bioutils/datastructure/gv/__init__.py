from __future__ import annotations

import uuid
from abc import abstractmethod
from typing import Callable

VALID_SORT_EXON_EXON_STRAND_POLICY = ("unstranded", "stranded", "none")
DEFAULT_SORT_EXON_EXON_STRAND_POLICY = "unstranded"
SequenceFuncType = Callable[[str, int, int], str]


class _NotSet:
    pass


_notset = _NotSet()


def generate_unknown_transcript_id() -> str:
    """Generate a new unknown transcript ID"""
    return 'unknown_transcript_id' + str(uuid.uuid4())


def generate_unknown_gene_id() -> str:
    """Generate a new unknown gene ID"""
    return 'unknown_gene_id' + str(uuid.uuid4())


class GVPError(ValueError):
    pass


class CanTranscribeInterface:
    @abstractmethod
    def transcribe(self, sequence_func: SequenceFuncType) -> str:
        raise NotImplementedError

    @abstractmethod
    @property
    def transcribed_length(self) -> int:
        raise NotImplementedError
