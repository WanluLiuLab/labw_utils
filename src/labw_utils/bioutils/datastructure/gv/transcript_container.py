from __future__ import annotations

from abc import abstractmethod
from typing import Iterable

from labw_utils.bioutils.datastructure.gv.transcript import Transcript


class TranscriptContainerInterface:

    @property
    @abstractmethod
    def number_of_transcripts(self) -> int:
        raise NotImplementedError

    @property
    @abstractmethod
    def transcript_values(self) -> Iterable[Transcript]:
        raise NotImplementedError

    @property
    @abstractmethod
    def transcript_ids(self) -> Iterable[str]:
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

