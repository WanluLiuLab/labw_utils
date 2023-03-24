from __future__ import annotations

__all__ = (
    "IntervalType",
    "NumpyIntervalEngine",
)

import functools
from collections import defaultdict
from typing import Iterable, Dict, Tuple, Optional
from typing import List

import numpy as np
import numpy.typing as npt
import pandas as pd

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_reader

IntervalType = Tuple[Tuple[str, Optional[bool]], int, int]


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


class NumpyIntervalEngine:
    """
    Store data in an NDArray with schema:

    [[s, e], [s, e], [s, e,], ...]
    """
    _chromosomal_split_np_index: Dict[Tuple[str, Optional[bool]], npt.NDArray]

    def _select_chromosome(
            self,
            query_chr: Tuple[str, Optional[bool]]
    ) -> Tuple[npt.NDArray, npt.NDArray]:
        stored_values_of_selected_chromosome = self._chromosomal_split_np_index[query_chr]
        s = stored_values_of_selected_chromosome[:, 0]
        e = stored_values_of_selected_chromosome[:, 1]
        return s, e

    def overlap(self, query_interval: IntervalType) -> Iterable[int]:
        query_chr, query_s, query_e = query_interval
        try:
            s, e = self._select_chromosome(query_chr)
        except KeyError:
            return None
        for it in np.nonzero(
                functools.reduce(
                    np.logical_or,
                    (
                            np.logical_and(
                                np.asarray(s < query_s),
                                np.asarray(query_s < e),
                            ),
                            np.logical_and(
                                np.asarray(s < query_e),
                                np.asarray(query_e < e),
                            ),
                            np.logical_and(
                                np.asarray(query_s < s),
                                np.asarray(s < query_e)
                            ),
                            np.logical_and(
                                np.asarray(query_s < e),
                                np.asarray(e < query_e)
                            )
                    )
                )
        )[0].tolist():
            yield it

    def __init__(self, chromosomal_split_np_index: Dict[Tuple[str, Optional[bool]], npt.NDArray]):
        self._chromosomal_split_np_index = chromosomal_split_np_index

    @classmethod
    def from_interval_tsv(cls, interval_file: str, show_tqdm: bool = True):
        pd_df = create_pandas_dataframe_from_input_file(
            interval_file, show_tqdm
        )
        tmpd = {}
        for chr_name in pd.unique(pd_df["chr"]):
            tmpd = pd_df.query(
                f"`chr` == '{chr_name}'"
            ).loc[:, ["s", "e"]].to_numpy(dtype=int)
        del pd_df
        return cls(tmpd)

    @classmethod
    def from_interval_iterator(cls, interval_iterator: Iterable[IntervalType]):
        tmpd: Dict[Tuple[str, Optional[bool]], List[Tuple[int, int]]] = defaultdict(lambda: [])
        for interval in interval_iterator:
            append_chr, append_s, append_e = interval
            tmpd[append_chr].append((append_s, append_e))
        return cls({k: np.array(tmpd[k], dtype=int) for k in tmpd.keys()})

    def match(self, query_interval: IntervalType) -> Iterable[int]:
        query_chr, query_s, query_e = query_interval
        try:
            s, e = self._select_chromosome(query_chr)
        except KeyError:
            return None
        match_result = np.nonzero(
            np.logical_and(
                np.asarray(s > query_s),
                np.asarray(e < query_e)
            )
        )[0]
        for it in match_result.tolist():
            yield it

    def __iter__(self) -> Iterable[IntervalType]:
        for chr_name, chr_value in self._chromosomal_split_np_index.items():
            for stored_values in chr_value:
                s, e = stored_values
                yield chr_name, s, e

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
