from collections import defaultdict
from typing import Iterable, Dict

import intervaltree

from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from naive_interval_engine import BaseNaiveIntervalEngine, IntervalType


class IntervalTreeIntervalEngine(BaseNaiveIntervalEngine):
    _dfs: Dict[str, intervaltree.IntervalTree]

    def overlap(self, query_interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = query_interval
        try:
            selected_chr = self._dfs[interval_chr]
        except KeyError:
            return None
        for it in sorted(map(
                lambda i: i.data,
                selected_chr.overlap(interval_s, interval_e)
        )):
            yield it

    def __init__(self, interval_file: str, show_tqdm: bool = True):
        if show_tqdm:
            reader = get_tqdm_line_reader(interval_file)
        else:
            reader = get_reader(interval_file)
        self._dfs = defaultdict(lambda: intervaltree.IntervalTree())
        _ = next(reader)
        for i, line in enumerate(reader):
            chr_name, s_str, e_str = line.rstrip("\n").split("\t")
            self._dfs[chr_name].addi(int(s_str), int(e_str), i)
        reader.close()

    def match(self, query_interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = query_interval
        try:
            selected_chr = self._dfs[interval_chr]
        except KeyError:
            return None
        for it in sorted(map(
                lambda i: i.data,
                selected_chr.envelop(interval_s, interval_e)
        )):
            yield it

    def __iter__(self) -> Iterable[IntervalType]:
        for chr_name, chr_value in self._dfs.items():
            for it in sorted(chr_value.items(), key=lambda it: it.data):
                yield chr_name, it.begin, it.end
