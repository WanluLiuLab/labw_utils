import numexpr as ne
import numpy as np

from labw_utils.typing_importer import Iterable
from naive_interval_engine import IntervalType
from naive_interval_engine.np_impl import NumpyIntervalEngine


class NumExprIntervalEngine(NumpyIntervalEngine):
    def overlap(self, query_interval: IntervalType) -> Iterable[int]:
        query_chr, query_s, query_e = query_interval
        try:
            s, e = self._select_chromosome(query_chr)
        except KeyError:
            return None
        for it in np.nonzero(
            ne.evaluate(
                "|".join(
                    (
                        "(" + "&".join(("(s < query_s)", "(query_s < e)")) + ")",
                        "(" + "&".join(("(s < query_e)", "(query_e < e)")) + ")",
                        "(" + "&".join(("(query_s < s)", "(s < query_e)")) + ")",
                        "(" + "&".join(("(query_s < e)", "(e < query_e)")) + ")",
                    )
                )
            )
        )[0].tolist():
            yield it

    def match(self, query_interval: IntervalType) -> Iterable[int]:
        query_chr, query_s, query_e = query_interval
        try:
            s, e = self._select_chromosome(query_chr)
        except KeyError:
            return None
        match_result = np.nonzero(ne.evaluate("(s > query_s) & (e < query_e)"))[0]
        for it in match_result.tolist():
            yield it
