from typing import Dict, Iterator, Iterable, Final

from labw_utils.bioutils.io import BaseFileIterator
from labw_utils.bioutils.typing.fai import FastaIndexEntry
from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader

FAI_INDEX_TYPE = Dict[str, FastaIndexEntry]


class FAIBasedFastaIndexIterator(BaseFileIterator):
    filetype: Final[str] = "FAI (FAI based)"

    def __init__(self, filename: str, show_tqdm: bool = True):
        super().__init__(filename, show_tqdm)
        if self._show_tqdm:
            self._fd = get_tqdm_line_reader(self.filename)
        else:
            self._fd = get_reader(self.filename)

    def __iter__(self) -> Iterator[FastaIndexEntry]:
        for line in self._fd:
            yield FastaIndexEntry.from_fai_str(line)
        self._fd.close()


class FastaBasedFastaIndexIterator(BaseFileIterator):
    filetype: Final[str] = "FAI (FASTA based)"
    _full_header: bool
    _name: str
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

    def _generate_fai_record(self) -> FastaIndexEntry:
        return FastaIndexEntry(
            name=self._name,
            length=self._length,
            offset=self._offset,
            line_blen=self._line_blen,
            line_len=self._line_len
        )

    def __iter__(self) -> Iterable[FastaIndexEntry]:
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
                    self._name = line[1:]
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


def create_fai_from_fasta(
        filename: str,
        index_filename: str
) -> FAI_INDEX_TYPE:
    """
    Create an FAI index for Fasta

    Do not use this feature on full headers.
    """
    return_fai: FAI_INDEX_TYPE = {}
    for fai_record in FastaBasedFastaIndexIterator(filename, full_header=False):
        return_fai[fai_record.name] = fai_record

    with get_writer(index_filename) as writer:
        for fai_record in return_fai.values():
            writer.write(str(fai_record) + "\n")
    return return_fai
