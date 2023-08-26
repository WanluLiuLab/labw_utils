"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

from abc import abstractmethod

from labw_utils.bioutils.datastructure.gv import GVPError
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.typing_importer import Sequence


class DuplicatedTranscriptIDError(GVPError):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    def __init__(self, transcript_id: str):
        super().__init__(f"Transcript ID {transcript_id} duplicated")


class TranscriptContainerInterface:
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    @property
    @abstractmethod
    def number_of_transcripts(self) -> int:
        raise NotImplementedError

    @property
    @abstractmethod
    def transcript_values(self) -> Sequence[Transcript]:
        raise NotImplementedError

    @property
    @abstractmethod
    def transcript_ids(self) -> Sequence[str]:
        raise NotImplementedError

    @abstractmethod
    def get_transcript(self, transcript_id: str) -> Transcript:
        raise NotImplementedError

    @abstractmethod
    def add_transcript(self, transcript: Transcript) -> TranscriptContainerInterface:
        raise NotImplementedError

    @abstractmethod
    def del_transcript(self, transcript_id: str) -> TranscriptContainerInterface:
        raise NotImplementedError

    @abstractmethod
    def replace_transcript(self, new_transcript: Transcript) -> TranscriptContainerInterface:
        raise NotImplementedError
