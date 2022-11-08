from __future__ import annotations

__all__ = (
    "parallel_map",
    "IntervalType",
    "BaseNaiveIntervalEngine"
)

from abc import abstractmethod
from typing import Tuple, List

IntervalType = Tuple[str, int, int]

import multiprocessing
from typing import Callable, TypeVar, Iterable

import joblib

_InType = TypeVar("_InType")
_OutType = TypeVar("_OutType")


def window(
        intervals: Iterable[_InType],
        batch_size: int,
) -> Iterable[List[_InType]]:
    cid = 0
    batch = []
    for interval in intervals:
        batch.append(interval)
        cid += 1
        if cid % batch_size == 0:
            yield batch
    if len(batch) != 0:
        yield batch


def parallel_map(
        f: Callable[[_InType], _OutType],
        input_iterable: Iterable[_InType],
        n_jobs: int = multiprocessing.cpu_count(),
        backend: str = "threading",
) -> Iterable[_OutType]:
    it: Iterable[_OutType] = joblib.Parallel(
        n_jobs=n_jobs,
        backend=backend
    )(
        joblib.delayed(f)(i) for i in input_iterable
    )
    return it


class BaseNaiveIntervalEngine:
    @abstractmethod
    def __init__(self, interval_file: str):
        _ = interval_file

    @abstractmethod
    def match(self, interval: IntervalType) -> Iterable[int]:
        pass

    def matches(self, intervals: Iterable[IntervalType]) -> Iterable[List[int]]:
        for interval in intervals:
            yield list(self.match(interval))

    def _parallel_prepared_matches(self, intervals: List[IntervalType]) -> List[List[int]]:
        return list(map(lambda x: list(self.match(x)), intervals))

    @abstractmethod
    def overlap(self, interval: IntervalType) -> Iterable[int]:
        pass

    def overlaps(self, intervals: Iterable[IntervalType]) -> Iterable[List[int]]:
        for interval in intervals:
            yield list(self.overlap(interval))

    def parallel_matches(
            self,
            intervals: Iterable[IntervalType],
            batch_size: int = 50
    ) -> Iterable[List[int]]:
        windows = list(window(intervals, batch_size))
        for matches in parallel_map(
                self._parallel_prepared_matches,
                windows,
                backend="loky"
        ):
            for matched in matches:
                yield matched

    def parallel_overlaps(
            self,
            intervals: Iterable[IntervalType],
            batch_size: int = 1000
    ) -> Iterable[List[int]]:
        for overlaps in parallel_map(
                self.overlaps,
                window(list(intervals), batch_size)
        ):
            for overlapped in overlaps:
                yield overlapped

    @abstractmethod
    def __iter__(self) -> Iterable[IntervalType]:
        pass
