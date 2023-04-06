from typing import Iterable, Iterator

from labw_utils.bioutils.parser import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.record.feature import Feature, DEFAULT_GTF_QUOTE_OPTIONS
from labw_utils.bioutils.record.gtf import parse_record, format_string
from labw_utils.commonutils.io.safe_io import get_writer
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

    def __init__(self, filename: str, quote: str = DEFAULT_GTF_QUOTE_OPTIONS, **kwargs):
        super().__init__(filename, **kwargs)
        self._fd = get_writer(self._filename)
        self._quote = quote

    @staticmethod
    def write_iterator(
            iterable: Iterable[Feature],
            filename: str,
            prefix_annotations: Iterable[str] = None,
            quote: str = DEFAULT_GTF_QUOTE_OPTIONS
    ):
        with GtfIteratorWriter(filename, quote) as writer:
            if prefix_annotations is not None:
                for annotation in prefix_annotations:
                    writer.write_comment(annotation)
            for feature in iterable:
                writer.write(feature)

    def write(self, record: Feature) -> None:
        self._fd.write(format_string(record, quote=self._quote) + "\n")

    def write_comment(self, comment: str):
        self._fd.write('#' + comment + "\n")
