import functools
from typing import Iterable, Dict

import numpy as np
import numpy.typing as npt
import pandas as pd

from naive_interval_engine import BaseNaiveIntervalEngine, IntervalType


class NumpyIntervalEngine(BaseNaiveIntervalEngine):
    _dfs: Dict[str, npt.NDArray]

    def overlap(self, interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = interval
        try:
            selected_chr = self._dfs[interval_chr]
        except KeyError:
            return None
        s = selected_chr[:, 0]
        e = selected_chr[:, 1]
        for it in np.nonzero(
                functools.reduce(
                    np.logical_or, (
                            np.logical_and(
                                np.asarray(s < interval_s),
                                np.asarray(interval_s < e),
                            ),
                            np.logical_and(
                                np.asarray(s < interval_e),
                                np.asarray(interval_e < e),
                            ),
                            np.logical_and(
                                np.asarray(interval_s < s),
                                np.asarray(s < interval_e)
                            ),
                            np.logical_and(
                                np.asarray(interval_s < e),
                                np.asarray(e < interval_e)
                            )
                    )
                )
        )[0].tolist():
            yield it

    def __init__(self, interval_file: str):
        pd_df = pd.read_csv(
            interval_file,
            sep="\t",
            engine="pyarrow",
            dtype={
                "chr": str,
                "s": int,
                "e": int
            }
        )
        self._dfs = {}
        for chr_name in pd.unique(pd_df["chr"]):
            self._dfs[chr_name] = pd_df.query(
                f"`chr` == '{chr_name}'"
            ).loc[:, ["s", "e"]].to_numpy(dtype=np.int_)

    def match(self, interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = interval
        try:
            selected_chr = self._dfs[interval_chr]
        except KeyError:
            return None
        s = selected_chr[:, 0]
        e = selected_chr[:, 1]
        match_result = np.nonzero(
            np.logical_and(
                np.asarray(s > interval_s),
                np.asarray(e < interval_e)
            )
        )[0]
        for it in match_result.tolist():
            yield it

    def __iter__(self) -> Iterable[IntervalType]:
        for chr_name, chr_value in self._dfs.items():
            for start_end in chr_value:
                yield chr_name, start_end[0], start_end[1]
