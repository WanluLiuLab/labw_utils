__all__ = (
    "FastaIndexParserError",
    "DuplicatedFastaNameError",
    "FastaBasedFastaIndexIterator",
    "FAIBasedFastaIndexIterator",
    "FastaIndexWriter"
)

from typing import Iterator, Iterable, Final, List

from labw_utils.bioutils.parser import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.record.fai import FastaIndexRecord
from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader


class FastaIndexParserError(ValueError):
    pass


class DuplicatedFastaNameError(FastaIndexParserError):
    def __init__(self, name: str):
        super().__init__(f"FAI seqname {name} had occurred more than once")


class FastaIndexNotWritableError(FastaIndexParserError):
    def __init__(self, name: str):
        super().__init__(
            f"FAI seqname '{name}' is valid in-memory not valid on disk"
            "Reason might be your using `full_header` on abnormal FASTAs"
        )


class FAIBasedFastaIndexIterator(BaseFileIterator):
    filetype: Final[str] = "FAI (FAI based)"

    def __init__(self, filename: str, show_tqdm: bool = True):
        super().__init__(filename, show_tqdm)
        if self._show_tqdm:
            self._fd = get_tqdm_line_reader(self.filename)
        else:
            self._fd = get_reader(self.filename)

    def __iter__(self) -> Iterator[FastaIndexRecord]:
        _names = []
        for line in self._fd:
            new_record = FastaIndexRecord.from_fai_str(line)
            if new_record.name in _names:
                raise DuplicatedFastaNameError(new_record.name)
            _names.append(new_record.name)
            yield new_record
        self._fd.close()


class FastaBasedFastaIndexIterator(BaseFileIterator):
    filetype: Final[str] = "FAI (FASTA based)"
    _full_header: bool
    _name: str
    _names: List[str]
    _length: int
    _offset: int
    _line_blen: int
    _line_len: int

    def __init__(
            self,
            filename: str,
            show_tqdm: bool = True,
            full_header: bool = True
    ):
        super().__init__(filename, show_tqdm)
        self._name = ""
        self._offset = 0
        self._line_len = 0
        self._line_blen = 0
        self._length = 0
        if self._show_tqdm:
            self._fd = get_tqdm_line_reader(self.filename)
        else:
            self._fd = get_reader(self.filename)
        self._full_header = full_header
        self._names = []

    def _generate_fai_record(self) -> FastaIndexRecord:
        if self._name in self._names:
            raise DuplicatedFastaNameError(self._name)
        else:
            self._names.append(self._name)
        return FastaIndexRecord(
            name=self._name,
            length=self._length,
            offset=self._offset,
            line_blen=self._line_blen,
            line_len=self._line_len
        )

    def __iter__(self) -> Iterable[FastaIndexRecord]:
        while True:
            line = self._fd.readline()
            if not line:
                break
            if line == '':
                continue
            if line[0] == '>':  # FASTA header
                if self._name != '':
                    yield self._generate_fai_record()
                    self._offset = 0
                    self._line_len = 0
                    self._line_blen = 0
                    self._length = 0
                if self._full_header:
                    self._name = line[1:].strip()
                else:
                    self._name = line[1:].strip().split(' ')[0].split('\t')[0]
                self._offset = self._fd.tell()
            else:
                if self._line_len == 0:
                    self._line_len = len(line)
                if self._line_blen == 0:
                    self._line_blen = len(line.rstrip('\r\n'))
                self._length += len(line.rstrip('\r\n'))
        if self._name != '':
            yield self._generate_fai_record()

        self._fd.close()


class FastaIndexWriter(BaseIteratorWriter):
    filetype: Final[str] = "FAI"

    def __init__(self, filename: str, **kwargs):
        super().__init__(filename, **kwargs)
        self._fd = get_writer(self._filename)

    def write(self, record: FastaIndexRecord) -> None:
        self._fd.write(str(record) + "\n")

    @staticmethod
    def write_iterator(iterable: Iterable[FastaIndexRecord], filename: str):
        with FastaIndexWriter(filename) as writer:
            for fastq_record in iterable:
                if fastq_record.name.count("\t") + fastq_record.name.count(" ") != 0:
                    writer.destroy_file()
                    raise FastaIndexNotWritableError(fastq_record.name)
                writer.write(fastq_record)
