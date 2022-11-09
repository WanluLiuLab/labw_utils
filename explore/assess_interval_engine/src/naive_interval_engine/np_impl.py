import functools
from typing import Iterable, Dict, Tuple

import numpy as np
import numpy.typing as npt
import pandas as pd

from naive_interval_engine import BaseNaiveIntervalEngine, IntervalType, create_pandas_dataframe_from_input_file


class NumpyIntervalEngine(BaseNaiveIntervalEngine):
    """
    Store data in an NDArray with schema:

    [[s, e], [s, e], [s, e,], ...]
    """
    _chromosomal_split_np_index: Dict[str, npt.NDArray[int]]

    def _select_chromosome(self, query_chr: str) -> Tuple[npt.NDArray[int], npt.NDArray[int]]:
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

    def __init__(self, interval_file: str, show_tqdm: bool = True):
        pd_df = create_pandas_dataframe_from_input_file(
            interval_file, show_tqdm
        )
        self._chromosomal_split_np_index = {}
        for chr_name in pd.unique(pd_df["chr"]):
            self._chromosomal_split_np_index[chr_name] = pd_df.query(
                f"`chr` == '{chr_name}'"
            ).loc[:, ["s", "e"]].to_numpy(dtype=np.int_)
        del pd_df

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
