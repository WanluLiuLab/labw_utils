"""
TODO: docs

.. versionadded:: 1.0.2
"""

from __future__ import annotations

from labw_utils import UnmetDependenciesError

__all__ = (
    "IntervalType",
    "NumpyIntervalEngine",
)

from collections import defaultdict
from labw_utils.typing_importer import Iterable, Dict, Tuple, Optional
from labw_utils.typing_importer import List

try:
    import numpy as np
    import numpy.typing as npt
except ImportError:
    raise UnmetDependenciesError("numpy")

IntervalType = Tuple[Tuple[str, Optional[bool]], int, int]


class NumpyIntervalEngine:
    """
    Store data in an NDArray with schema:

    [[s, e], [s, e], [s, e,], ...]

    .. versionadded:: 1.0.2
    """

    _chromosomal_split_np_index: Dict[Tuple[str, Optional[bool]], npt.NDArray]
    _interval_name_index: Dict[IntervalType, str]

    def _select_chromosome(self, query_chr: Tuple[str, Optional[bool]]) -> Tuple[npt.NDArray, npt.NDArray]:
        stored_values_of_selected_chromosome = self._chromosomal_split_np_index[query_chr]
        s = stored_values_of_selected_chromosome[:, 0]
        e = stored_values_of_selected_chromosome[:, 1]
        return s, e

    def overlap(self, query_interval: IntervalType) -> Iterable[str]:
        query_chr, query_s, query_e = query_interval
        try:
            s, e = self._select_chromosome(query_chr)
        except KeyError:
            return None
        es = np.ndarray((2, e.shape[0]))
        es[0:] = query_e
        es[1:] = e
        ss = np.ndarray((2, s.shape[0]))
        ss[0:] = query_s
        ss[1:] = s
        for match_id in np.where(np.min(es, axis=0) - np.max(ss, axis=0) > 0)[0].tolist():
            yield self._interval_name_index[(query_chr, s[match_id], e[match_id])]

    def __init__(
        self,
        chromosomal_split_np_index: Dict[Tuple[str, Optional[bool]], npt.NDArray],
        interval_name_index: Dict[IntervalType, str],
    ):
        self._chromosomal_split_np_index = chromosomal_split_np_index
        self._interval_name_index = interval_name_index

    @classmethod
    def from_interval_iterator(cls, interval_iterator: Iterable[IntervalType], names_iterator: Iterable[str]):
        tmpd: Dict[Tuple[str, Optional[bool]], List[Tuple[int, int]]] = defaultdict(lambda: [])
        interval_name_index: Dict[IntervalType, str] = {}
        for interval, name in zip(interval_iterator, names_iterator):
            append_chr, append_s, append_e = interval
            tmpd[append_chr].append((append_s, append_e))
            interval_name_index[interval] = name
        return cls({k: np.array(tmpd[k], dtype=int) for k in tmpd.keys()}, interval_name_index)
