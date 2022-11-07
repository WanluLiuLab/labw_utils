from typing import Iterator, Union, Final

from labw_utils.bioutils.io import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.typing.fastq import FastqRecord
from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader


class FastqIterator(BaseFileIterator):
    filetype: Final[str] = "FASTQ"

    def __init__(self, filename: str, show_tqdm: bool = True):
        super().__init__(filename, show_tqdm)
        if self._show_tqdm:
            self._fd = get_tqdm_line_reader(self.filename)
        else:
            self._fd = get_reader(self.filename)

    def __iter__(self) -> Iterator[FastqRecord]:
        while True:
            lines = [self._fd.readline(-1) for _ in range(4)]
            if '' in lines:
                break
            yield FastqRecord.from_str(lines)
        self._fd.close()


class FastqWriter(BaseIteratorWriter):
    filetype: Final[str] = "FASTQ"

    @staticmethod
    def write_iterator(
            iterable: Union[FastqIterator, Iterator[FastqRecord]],
            filename: str
    ):
        with FastqWriter(filename) as writer:
            for fastq_record in iterable:
                writer.write(fastq_record)

    def __init__(self, filename: str, **kwargs):
        super().__init__(filename, **kwargs)
        self._fd = get_writer(self._filename)

    def write(self, fastq_record: FastqRecord):
        self._fd.write(str(fastq_record) + "\n")
