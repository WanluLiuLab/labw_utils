from typing import Iterable, Dict

import numexpr as ne
import numpy as np
import numpy.typing as npt
import pandas as pd

from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_reader
from naive_interval_engine import BaseNaiveIntervalEngine, IntervalType


class NumExprIntervalEngine(BaseNaiveIntervalEngine):
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
                ne.evaluate(
                    "|".join((
                            "(" + "&".join((
                                    "(s < interval_s)",
                                    "(interval_s < e)"
                            )) + ")",
                            "(" + "&".join((
                                    "(s < interval_e)",
                                    "(interval_e < e)"
                            )) + ")",
                            "(" + "&".join((
                                    "(interval_s < s)",
                                    "(s < interval_e)"
                            )) + ")",
                            "(" + "&".join((
                                    "(interval_s < e)",
                                    "(e < interval_e)"
                            )) + ")",
                    ))
                )
        )[0].tolist():
            yield it

    def __init__(self, interval_file: str, show_tqdm: bool = True):
        if show_tqdm:
            reader = get_tqdm_reader(interval_file, is_binary=True)
        else:
            reader = get_reader(interval_file, is_binary=True)
        pd_df = pd.read_csv(
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
        self._dfs = {}
        for chr_name in pd.unique(pd_df["chr"]):
            self._dfs[chr_name] = pd_df.query(
                f"`chr` == '{chr_name}'"
            ).loc[:, ["s", "e"]].to_numpy(dtype=np.int_)
        del pd_df

    def match(self, interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = interval
        try:
            selected_chr = self._dfs[interval_chr]
        except KeyError:
            return None
        s = selected_chr[:, 0]
        e = selected_chr[:, 1]
        match_result = np.nonzero(
            ne.evaluate(
                "(s > interval_s) & (e < interval_e)"
            )
        )[0]
        for it in match_result.tolist():
            yield it

    def __iter__(self) -> Iterable[IntervalType]:
        for chr_name, chr_value in self._dfs.items():
            for start_end in chr_value:
                yield chr_name, start_end[0], start_end[1]
