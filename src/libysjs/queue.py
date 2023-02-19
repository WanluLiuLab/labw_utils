from __future__ import annotations

from typing import List, Mapping, Union, Optional


class YSJSDQueueStatus:
    ...


class YSJSDQueueConfig:
    _queue_name: str
    _queue_description: str
    _queue_max_capacity: int

    def __init__(
            self,
            queue_name: str,
            queue_description: str,
            queue_max_capacity: int
    ):
        self._queue_name = queue_name
        self._queue_description = queue_description
        self._queue_max_capacity = queue_max_capacity

    def to_dict(self) -> Mapping[str, Union[str, int]]:
        return {
            "queue_name": self._queue_name,
            "queue_description": self._queue_description,
            "queue_max_capacity": self._queue_max_capacity
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Union[str, int]]):
        return cls(**in_dict)

