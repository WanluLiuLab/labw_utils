from __future__ import annotations

import uuid
from abc import abstractmethod, ABC

from labw_utils.typing_importer import Callable

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


class CanTranscribeInterface(ABC):
    __slots__ = []

    @abstractmethod
    def transcribe(self, sequence_func: SequenceFuncType) -> str:
        raise NotImplementedError

    @property
    @abstractmethod
    def transcribed_length(self) -> int:
        raise NotImplementedError


class SortedContainerInterface:
    _is_sorted: bool

    @property
    def is_sorted(self) -> bool:
        return self._is_sorted


class CanCheckInterface:
    _is_checked: bool

    @property
    def is_checked(self) -> bool:
        return self._is_checked
