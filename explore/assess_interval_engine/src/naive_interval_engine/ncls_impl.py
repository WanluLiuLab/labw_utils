"""
Have bugs; do not use.

See: <https://github.com/biocore-ntnu/ncls/issues/43>

Also waiting for more documentations.
"""

from enum import Enum

import ncls
import pandas as pd

from labw_utils.typing_importer import Iterable, Dict, List, Optional
from naive_interval_engine import IntervalEngineType, IntervalType, create_pandas_dataframe_from_input_file


class NclsIntervalEngine(IntervalEngineType):
    _chromosomal_split_ncls_index: Dict[str, ncls.NCLS64]

    class _OperationType(Enum):
        Match = 0
        Overlap = 1

    def __init__(self, interval_file: str, show_tqdm: bool = True):
        self._chromosomal_split_ncls_index = {}
        pd_df = create_pandas_dataframe_from_input_file(
            interval_file, show_tqdm
        )
        for chr_name in pd.unique(pd_df["chr"]):
            pd_df_of_selected_chromosome = self._chromosomal_split_ncls_index[chr_name] = pd_df.query(
                f"`chr` == '{chr_name}'"
            ).loc[:, ["s", "e", "idx"]]
            self._chromosomal_split_ncls_index[chr_name] = ncls.NCLS(
                starts=pd_df_of_selected_chromosome["s"].values,
                ends=pd_df_of_selected_chromosome["e"].values,
                ids=pd_df_of_selected_chromosome["idx"].values,
            )
        del pd_df

    def _matches_or_overlaps(
            self,
            operation_type: _OperationType,
            query_intervals: Iterable[IntervalType],
            show_tqdm: bool = True
    ) -> Iterable[List[int]]:
        _ = show_tqdm

        query_intervals_df = pd.DataFrame(data=list(query_intervals), columns=["chr", "s", "e"])
        """
        Query dataframe. Schema:

        - chr, s, e: same.
        - qidx: Index used for sorting after querying.
        """

        full_result_df: Optional[pd.DataFrame] = None
        """
        Result dataframe. Schema:

        - qidx
        - idx
        """

        query_length = query_intervals_df.shape[0]
        query_intervals_df['qidx'] = range(0, query_length)

        if operation_type == self._OperationType.Match:
            operation_func = ncls.NCLS64.all_containments_both
        else:
            operation_func = ncls.NCLS64.all_overlaps_both

        for chr_name in pd.unique(query_intervals_df["chr"]):
            if chr_name not in self._chromosomal_split_ncls_index:
                continue
            query_df = query_intervals_df.query(
                f"`chr` == '{chr_name}'"
            ).loc[:, ["s", "e", "qidx"]]
            result_qidx, result_idx = operation_func(
                self._chromosomal_split_ncls_index[chr_name],
                query_df["s"].values,
                query_df["e"].values,
                query_df["qidx"].values,
            )
            result_df_of_selected_chromosome = pd.DataFrame({
                "qidx": result_qidx,
                "idx": result_idx
            })
            if full_result_df is None:
                full_result_df = result_df_of_selected_chromosome
            else:
                full_result_df = pd.concat([full_result_df, result_df_of_selected_chromosome])
        if full_result_df is None:
            full_result_df = pd.DataFrame({
                "qidx": [],
                "idx": []
            })
        for qidx in range(query_length):
            result_df_of_selected_qidx = full_result_df.query(
                f"`qidx` == {qidx}"
            )
            yield list(sorted(result_df_of_selected_qidx["idx"]))

    def overlap(self, query_interval: IntervalType) -> Iterable[int]:
        return list(self.overlaps(
            [query_interval],
        ))[0]

    def match(self, query_interval: IntervalType) -> Iterable[int]:
        return list(self.matches(
            [query_interval],
        ))[0]

    def matches(self, query_intervals: Iterable[IntervalType], show_tqdm: bool = True) -> Iterable[List[int]]:
        return self._matches_or_overlaps(self._OperationType.Match, query_intervals=query_intervals,
                                         show_tqdm=show_tqdm)

    def overlaps(self, query_intervals: Iterable[IntervalType], show_tqdm: bool = True) -> Iterable[List[int]]:
        return self._matches_or_overlaps(self._OperationType.Overlap, query_intervals=query_intervals,
                                         show_tqdm=show_tqdm)

    def __iter__(self) -> Iterable[IntervalType]:
        for chr_name, ncls_index_of_selected_chromosome in self._chromosomal_split_ncls_index.items():
            for it in sorted(
                    ncls_index_of_selected_chromosome.intervals(),
                    key=lambda it: it[2]
            ):
                s, e, _ = it
                yield chr_name, s, e
