import pandas as pd

from labw_utils.typing_importer import Iterable
from naive_interval_engine import BaseNaiveIntervalEngine, IntervalType, create_pandas_dataframe_from_input_file


class PandasIntervalEngine(BaseNaiveIntervalEngine):
    _df: pd.DataFrame

    def overlap(self, query_interval: IntervalType) -> Iterable[int]:
        query_chr, query_s, query_e = query_interval
        for idx in self._df.query(
                f"`chr` == '{query_chr}' & (" +
                "|".join((
                        f"s < {query_s} < e",
                        f"s < {query_e} < e",
                        f"{query_s} < s < {query_e}",
                        f"{query_s} < e < {query_e}"
                )) +
                ")"
        )["idx"]:
            yield idx

    def __init__(self, interval_file: str, show_tqdm: bool = True):
        self._df = create_pandas_dataframe_from_input_file(
            interval_file, show_tqdm
        )

    def match(self, query_interval: IntervalType) -> Iterable[int]:
        query_chr, query_s, query_e = query_interval
        for idx in self._df.query(
                f"`chr` == '{query_chr}' & s > {query_s} & e < {query_e}"
        )["idx"]:
            yield idx

    def __iter__(self) -> Iterable[IntervalType]:
        for it in self._df.itertuples(index=False):
            yield it[0:3]
