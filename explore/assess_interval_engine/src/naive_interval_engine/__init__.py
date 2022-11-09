from __future__ import annotations

__all__ = (
    "IntervalType",
    "IntervalEngineType",
    "BaseNaiveIntervalEngine",
)

from abc import abstractmethod, ABC
from typing import Tuple, List, Iterable

import pandas as pd

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_reader

IntervalType = Tuple[str, int, int]


def create_pandas_dataframe_from_input_file(
        interval_file: str,
        show_tqdm: bool = True
) -> pd.DataFrame:
    """
    Pandas dataframe for the intervals.

    Schema:
    - chr: str, chromosome name
    - s: int, sequence start
    - e: int, sequence end
    - idx: int, an auto-incremental index
    """
    if show_tqdm:
        reader = get_tqdm_reader(interval_file, is_binary=True)
    else:
        reader = get_reader(interval_file, is_binary=True)
    df = pd.read_csv(
        reader,
        sep="\t",
        engine="pyarrow",
        dtype={
            "chr": str,
            "s": int,
            "e": int
        }
    )
    reader.close()
    df["idx"] = range(0, df.shape[0])
    return df


class IntervalEngineType:
    @abstractmethod
    def __init__(
            self,
            interval_file: str,
            show_tqdm: bool = True
    ):
        _ = interval_file, show_tqdm

    @abstractmethod
    def match(self, query_interval: IntervalType) -> Iterable[int]:
        pass

    @abstractmethod
    def matches(
            self,
            query_intervals: Iterable[IntervalType],
            show_tqdm: bool = True
    ) -> Iterable[List[int]]:
        pass

    @abstractmethod
    def overlap(self, query_interval: IntervalType) -> Iterable[int]:
        pass

    @abstractmethod
    def overlaps(
            self,
            query_intervals: Iterable[IntervalType],
            show_tqdm: bool = True
    ) -> Iterable[List[int]]:
        pass

    @abstractmethod
    def __iter__(self) -> Iterable[IntervalType]:
        pass


class BaseNaiveIntervalEngine(IntervalEngineType, ABC):
    def matches(self, query_intervals: Iterable[IntervalType], show_tqdm: bool = True) -> Iterable[List[int]]:
        if show_tqdm:
            query_intervals = tqdm(
                iterable=list(query_intervals),
                desc="matching...",
            )
        for interval in query_intervals:
            yield list(self.match(interval))

    def overlaps(self, query_intervals: Iterable[IntervalType], show_tqdm: bool = True) -> Iterable[List[int]]:
        if show_tqdm:
            query_intervals = tqdm(
                iterable=list(query_intervals),
                desc="overlapping...",
            )
        for interval in query_intervals:
            yield list(self.overlap(interval))
