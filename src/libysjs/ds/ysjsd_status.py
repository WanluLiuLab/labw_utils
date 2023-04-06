from typing import Any, Mapping

AVAILABLE_YSJSD_STATE = ("starting", "running", "terminating")


class YSJSDStatus:
    _state: str
    _current_cpu: float
    _current_mem: float
    _pending_queue_length: int
    _running_queue_length: int
    _uptime: float

    def __init__(
            self,
            state: str,
            current_cpu: float,
            current_mem: float,
            pending_queue_length: int,
            running_queue_length: int,
            uptime: float
    ) -> None:
        self._state = state
        self._current_cpu = current_cpu
        self._current_mem = current_mem
        self._pending_queue_length = pending_queue_length
        self._running_queue_length = running_queue_length
        self._uptime = uptime

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "state": self._state,
            "current_cpu": self._current_cpu,
            "current_mem": self._current_mem,
            "pending_queue_length": self._pending_queue_length,
            "running_queue_length": self._running_queue_length,
            "uptime": self._uptime
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

    @property
    def state(self) -> str:
        return self._state

    @property
    def current_cpu(self) -> float:
        return self._current_cpu

    @property
    def current_mem(self) -> float:
        return self._current_mem

    @property
    def pending_queue_length(self) -> int:
        return self._pending_queue_length

    @property
    def running_queue_length(self) -> int:
        return self._running_queue_length

    @property
    def uptime(self) -> float:
        return self._uptime
