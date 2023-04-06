__all__ = (
    "FastaIterator",
    "extract_fasta_name",
    "FastaWriter"
)

from typing import Iterable, Final

from labw_utils.bioutils.parser import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader


def extract_fasta_name(line: str, full_header: bool) -> str:
    if full_header:
        return line[1:].strip()
    else:
        return line[1:].strip().split(' ')[0].split('\t')[0]


class FastaIterator(BaseFileIterator):
    filetype: Final[str] = "FASTA"
    _full_header: bool

    def __init__(
            self,
            filename: str,
            show_tqdm: bool = True,
            full_header: bool = True
    ):
        super().__init__(filename, show_tqdm)
        if self._show_tqdm:
            # Here required binary to deal with line endings
            self._fd = get_tqdm_line_reader(self.filename, is_binary=True)
        else:
            self._fd = get_reader(self.filename, is_binary=True)
        self._full_header = full_header

    def __iter__(self) -> Iterable[FastaRecord]:
        chr_name = ""
        seq = ""
        if self._show_tqdm:
            it = get_tqdm_line_reader(self.filename)
        else:
            it = get_reader(self.filename)
        for line in it:
            if line == "":
                continue
            if line[0] == '>':  # FASTA header
                if chr_name != '':
                    yield FastaRecord(
                        seq_id=chr_name,
                        sequence=seq
                    )
                    seq = ''
                if self._full_header:
                    chr_name = line[1:].strip()
                else:
                    chr_name = line[1:].strip().split(' ')[0].split('\t')[0]
            else:
                seq = seq + line.strip()
        if chr_name != '':
            yield FastaRecord(
                seq_id=chr_name,
                sequence=seq
            )
        self._fd.close()


class FastaWriter(BaseIteratorWriter):
    filetype: Final[str] = "FASTA"

    def __init__(self, filename: str, **kwargs):
        super().__init__(filename, **kwargs)
        self._fd = get_writer(self._filename)

    def write(self, record: FastaRecord) -> None:
        self._fd.write(str(record) + "\n")

    @staticmethod
    def write_iterator(iterable: Iterable[FastaRecord], filename: str):
        with FastaWriter(filename) as writer:
            for fastq_record in iterable:
                writer.write(fastq_record)
