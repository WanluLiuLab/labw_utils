from __future__ import annotations

from labw_utils.typing_importer import Optional, Mapping, Any
from libysjs.ds.ysjs_submission import YSJSSubmission

AVAILABLE_JOB_STATUS = (
    "pending",
    "starting",
    "running",
    "canceled",
    "finished",
)


class YSJSJob:
    _submission: YSJSSubmission
    _status: str
    _retv: Optional[int]
    _job_id: int
    _start_time: Optional[float]
    _terminate_time: Optional[float]
    _pid: Optional[int]

    def __init__(
        self,
        submission: YSJSSubmission,
        job_id: int,
        status: str,
        retv: Optional[int],
        start_time: Optional[float],
        terminate_time: Optional[float],
        pid: Optional[int],
    ):
        self._submission = submission
        self._status = status
        self._retv = retv
        self._job_id = job_id
        self._start_time = start_time
        self._terminate_time = terminate_time
        self._pid = pid

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "submission": self._submission,
            "status": self._status,
            "retv": self._retv,
            "job_id": self._job_id,
            "start_time": self._start_time,
            "terminate_time": self._terminate_time,
            "pid": self._pid,
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

    @classmethod
    def new(cls, submission: YSJSSubmission, job_id: int, **kwargs):
        return cls(
            submission=submission,
            status="pending",
            retv=None,
            job_id=job_id,
            start_time=None,
            terminate_time=None,
            pid=None,
        )

    @property
    def submission(self) -> YSJSSubmission:
        return self._submission

    @property
    def job_id(self) -> int:
        return self._job_id

    @property
    def status(self) -> str:
        return self._status

    @property
    def start_time(self) -> Optional[float]:
        return self._start_time

    def __repr__(self):
        return f"job (id={self._job_id}, status={self._status})"
