from typing import Iterable

import pandas as pd

from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_reader
from naive_interval_engine import BaseNaiveIntervalEngine, IntervalType


class PandasIntervalEngine(BaseNaiveIntervalEngine):
    _pd: pd.DataFrame

    def overlap(self, interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = interval
        for idx in self._pd.query(
                f"`chr` == '{interval_chr}' & (" +
                "|".join((
                        f"s < {interval_s} < e",
                        f"s < {interval_e} < e",
                        f"{interval_s} < s < {interval_e}",
                        f"{interval_s} < e < {interval_e}"
                )) +
                ")"
        )["idx"]:
            yield idx

    def __init__(self, interval_file: str, show_tqdm: bool = True):
        if show_tqdm:
            reader = get_tqdm_reader(interval_file, is_binary=True)
        else:
            reader = get_reader(interval_file, is_binary=True)
        self._pd = pd.read_csv(
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
        self._pd["idx"] = range(0, self._pd.shape[0])

    def match(self, interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = interval
        for idx in self._pd.query(
                f"`chr` == '{interval_chr}' & s > {interval_s} & e < {interval_e}"
        )["idx"]:
            yield idx

    def __iter__(self) -> Iterable[IntervalType]:
        for it in self._pd.itertuples(index=False):
            yield it[0:3]
