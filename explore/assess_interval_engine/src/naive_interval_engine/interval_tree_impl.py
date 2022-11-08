from collections import defaultdict
from typing import Iterable, Dict

import intervaltree

from naive_interval_engine import BaseNaiveIntervalEngine, IntervalType


class IntervalTreeIntervalEngine(BaseNaiveIntervalEngine):
    _dfs: Dict[str, intervaltree.IntervalTree]

    def overlap(self, interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = interval
        try:
            selected_chr = self._dfs[interval_chr]
        except KeyError:
            return None
        for it in sorted(map(
                lambda i: i.data,
                selected_chr.overlap(interval_s, interval_e)
        )):
            yield it

    def __init__(self, interval_file: str):
        self._dfs = defaultdict(lambda: intervaltree.IntervalTree())
        with open(interval_file, "rt") as reader:
            _ = next(reader)
            for i, line in enumerate(reader):
                chr_name, s_str, e_str = line.rstrip("\n").split("\t")
                self._dfs[chr_name].addi(int(s_str), int(e_str), i)

    def match(self, interval: IntervalType) -> Iterable[int]:
        interval_chr, interval_s, interval_e = interval
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
