from __future__ import annotations

import bisect
from abc import abstractmethod
from typing import List, Iterable

from labw_utils.bioutils.datastructure.gv import GVPError
from labw_utils.bioutils.datastructure.gv.transcript import Transcript


class DuplicatedTranscriptIDError(GVPError):
    pass


class TranscriptContainerInterface:

    @abstractmethod
    @property
    def number_of_transcripts(self) -> int:
        raise NotImplementedError

    @abstractmethod
    @property
    def transcript_values(self) -> Iterable[Transcript]:
        raise NotImplementedError

    @abstractmethod
    @property
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


class TranscriptContainer(TranscriptContainerInterface):
    __slots__ = (
        "_transcripts",
        "_transcript_ids"
    )
    _transcripts: List[Transcript]
    _transcript_ids: List[str]

    @property
    def number_of_transcripts(self) -> int:
        return len(self._transcripts)

    @property
    def transcript_values(self) -> Iterable[Transcript]:
        return iter(self._transcripts)

    @property
    def transcript_ids(self) -> Iterable[str]:
        return iter(self._transcript_ids)

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self._transcripts[self._transcript_ids.index(transcript_id)]

    def __init__(
            self,
            transcripts: Iterable[Transcript]
    ):
        self._transcripts = []
        self._transcript_ids = []
        self._transcripts = list(transcripts)
        self._transcript_ids = list(transcript.transcript_id for transcript in transcripts)

    def add_transcript(self, transcript: Transcript) -> TranscriptContainer:
        new_transcripts = list(self._transcripts)
        if transcript.transcript_id in self._transcripts:
            raise DuplicatedTranscriptIDError(
                f"Transcript ID {transcript.transcript_id} duplicated"
            )
        new_pos = bisect.bisect_left(self._transcripts, transcript)
        new_transcripts.insert(new_pos, transcript)
        return TranscriptContainer(new_transcripts)

    def del_transcript(self, transcript_id: str) -> TranscriptContainer:
        new_transcripts = list(self._transcripts)
        _ = new_transcripts.pop(self._transcript_ids.index(transcript_id))
        return TranscriptContainer(new_transcripts)

    def replace_transcript(self, new_transcript: Transcript) -> TranscriptContainer:
        return self.del_transcript(new_transcript.transcript_id).add_transcript(new_transcript)
