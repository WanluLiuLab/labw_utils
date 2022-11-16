from typing import Iterable, Iterator

from labw_utils.bioutils.parser import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.record.feature import Feature
from labw_utils.bioutils.record.gtf import parse_record, format_string
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader


class GtfIterator(BaseFileIterator, Iterable[Feature]):
    filetype: str = "GTF"
    record_type = Feature

    def __iter__(self) -> Iterator[Feature]:
        for line in get_tqdm_line_reader(self.filename):
            if line.startswith('#') or line == '':
                continue
            yield parse_record(line)


class GtfIteratorWriter(BaseIteratorWriter):
    filetype: str = "GTF"
    record_type = Feature

    @staticmethod
    def write_iterator(
            iterable: Iterable[Feature],
            filename: str,
            prefix_annotations: Iterable[str] = None
    ):
        with GtfIteratorWriter(filename) as writer:
            if prefix_annotations is not None:
                for annotation in prefix_annotations:
                    writer.write_comment(annotation)
            for feature in iterable:
                writer.write_feature(feature)

    def write(self, record: Feature) -> None:
        self._fd.write(format_string(record) + "\n")

    def write_comment(self, comment: str):
        self._fd.write('#' + comment + "\n")
