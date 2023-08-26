"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

import uuid
from abc import abstractmethod, ABC

from labw_utils.typing_importer import Callable, List, Tuple

VALID_SORT_EXON_EXON_STRAND_POLICY = ("unstranded", "stranded", "none")
"""
TODO: docs

.. versionadded:: 1.0.2
"""

DEFAULT_SORT_EXON_EXON_STRAND_POLICY = "unstranded"
"""
TODO: docs

.. versionadded:: 1.0.2
"""

SequenceFuncType = Callable[[str, int, int], str]
"""
TODO: docs

.. versionadded:: 1.0.2
"""

LegalizeRegionFuncType = Callable[[str, int, int], Tuple[str, int, int]]
"""
TODO: docs

.. versionadded:: 1.0.2
"""


def dumb_legalize_region_func(seqname: str, from_pos: int, to_pos: int) -> Tuple[str, int, int]:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    return seqname, from_pos, to_pos


def generate_unknown_transcript_id() -> str:
    """
    Generate a new unknown transcript ID

    .. versionadded:: 1.0.2
    """
    return "unknown_transcript_id" + str(uuid.uuid4())


def generate_unknown_gene_id() -> str:
    """
    Generate a new unknown gene ID

    .. versionadded:: 1.0.2
    """
    return "unknown_gene_id" + str(uuid.uuid4())


class GVPError(ValueError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    pass


class CanTranscribeInterface(ABC):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    __slots__: List[str] = []

    @abstractmethod
    def transcribe(
        self, sequence_func: SequenceFuncType, legalize_region_func: LegalizeRegionFuncType = dumb_legalize_region_func
    ) -> str:
        raise NotImplementedError

    @property
    @abstractmethod
    def transcribed_length(self) -> int:
        raise NotImplementedError


class CanCheckInterface:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    _is_checked: bool

    @property
    def is_checked(self) -> bool:
        return self._is_checked
