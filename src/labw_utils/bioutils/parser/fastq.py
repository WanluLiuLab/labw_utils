"""
TODO: docs

.. versionadded:: 1.0.2
"""

__all__ = ("FastqIterator", "FastqWriter")

from labw_utils.bioutils.parser import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.record.fastq import FastqRecord
from labw_utils.commonutils.lwio.safe_io import get_writer, get_reader
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.typing_importer import Iterator, Final, Iterable


class FastqIterator(BaseFileIterator, Iterable[FastqRecord]):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    filetype: Final[str] = "FASTQ"
    _full_header: bool

    def __init__(self, filename: str, show_tqdm: bool = True, full_header: bool = True):
        super().__init__(filename, show_tqdm)
        if self._show_tqdm:
            self._fd = get_tqdm_line_reader(self.filename)
        else:
            self._fd = get_reader(self.filename)
        self._full_header = full_header

    def __iter__(self) -> Iterator[FastqRecord]:
        while True:
            lines = [self._fd.readline(-1) for _ in range(4)]
            if "" in lines or len(lines) != 4:
                break
            yield FastqRecord.from_str(lines, full_header=self._full_header)
        self._fd.close()


class FastqWriter(BaseIteratorWriter):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    filetype: Final[str] = "FASTQ"

    @staticmethod
    def write_iterator(iterable: Iterator[FastqRecord], filename: str, **kwargs):
        with FastqWriter(filename) as writer:
            for fastq_record in iterable:
                writer.write(fastq_record)

    def __init__(self, filename: str, **kwargs):
        super().__init__(filename, **kwargs)
        self._fd = get_writer(self._filename)

    def write(self, record: FastqRecord):
        self._fd.write(str(record) + "\n")
