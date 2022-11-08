from __future__ import annotations

__all__ = (
    "IntervalType",
    "IntervalEngineType",
    "BaseNaiveIntervalEngine",
)

from abc import abstractmethod, ABC
from typing import Tuple, List, Iterable

from labw_utils.commonutils.importer.tqdm_importer import tqdm

IntervalType = Tuple[str, int, int]


class IntervalEngineType:
    @abstractmethod
    def __init__(
            self,
            interval_file: str,
            show_tqdm: bool = True
    ):
        _ = interval_file, show_tqdm

    @abstractmethod
    def match(self, interval: IntervalType) -> Iterable[int]:
        pass

    @abstractmethod
    def matches(
            self,
            intervals: Iterable[IntervalType],
            show_tqdm: bool = True
    ) -> Iterable[List[int]]:
        pass

    @abstractmethod
    def overlap(self, interval: IntervalType) -> Iterable[int]:
        pass

    @abstractmethod
    def overlaps(
            self,
            intervals: Iterable[IntervalType],
            show_tqdm: bool = True
    ) -> Iterable[List[int]]:
        pass

    @abstractmethod
    def __iter__(self) -> Iterable[IntervalType]:
        pass


class BaseNaiveIntervalEngine(IntervalEngineType, ABC):
    def matches(
            self,
            intervals: Iterable[IntervalType],
            show_tqdm: bool = True
    ) -> Iterable[List[int]]:
        if show_tqdm:
            intervals = tqdm(
                iterable=list(intervals),
                desc="matching...",
            )
        for interval in intervals:
            yield list(self.match(interval))

    def overlaps(
            self,
            intervals: Iterable[IntervalType],
            show_tqdm: bool = True
    ) -> Iterable[List[int]]:
        if show_tqdm:
            intervals = tqdm(
                iterable=list(intervals),
                desc="overlapping...",
            )
        for interval in intervals:
            yield list(self.overlap(interval))
