from __future__ import annotations

from typing import Optional

from libysjs.submission import YSJSSubmission


class CoreYSJSJob:
    _submission: YSJSSubmission
    _status: str
    _retv: Optional[int]

    def __init__(self, submission: YSJSSubmission):
        super().__init__()
        self._submission = submission
        self._status = "pending"
        self._retv = None

    @property
    def submission(self) -> YSJSSubmission:
        return self._submission

    @property
    def status(self) -> str:
        return self._status
