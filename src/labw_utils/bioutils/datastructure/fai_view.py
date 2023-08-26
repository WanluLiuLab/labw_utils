"""
TODO: docs

 .. versionadded:: 1.0.2
"""
from __future__ import annotations

__all__ = ("FastaIndexView", "create_fai_from_fasta")

from labw_utils.bioutils.parser.fai import FAIBasedFastaIndexIterator, FastaBasedFastaIndexIterator, FastaIndexWriter
from labw_utils.bioutils.record.fai import FastaIndexRecord
from labw_utils.typing_importer import Dict


class FastaIndexView:
    """
    TODO: docs

     .. versionadded:: 1.0.2
    """

    _d: Dict[str, FastaIndexRecord]
    _full_header: bool
    _filename: str

    def __init__(self, filename: str, d: Dict[str, FastaIndexRecord], full_header: bool):
        self._filename = filename
        self._d = d
        self._full_header = full_header

    @classmethod
    def from_fai(cls, filename: str, show_tqdm: bool = True):
        with FAIBasedFastaIndexIterator(filename, show_tqdm=show_tqdm) as faii:
            d = {fai_record.name: fai_record for fai_record in faii}
            return cls(filename, d, False)

    @classmethod
    def from_fasta(cls, filename: str, full_header: bool = True, show_tqdm: bool = True):
        with FastaBasedFastaIndexIterator(filename, show_tqdm=show_tqdm, full_header=full_header) as faii:
            d = {fai_record.name: fai_record for fai_record in faii}
            return cls(filename, d, True)

    def write(self, filename: str) -> None:
        FastaIndexWriter.write_iterator(self._d.values(), filename)

    def keys(self):
        return self._d.keys()

    def values(self):
        return self._d.values()

    def __getitem__(self, item):
        return self._d[item]

    def __len__(self):
        return len(self._d)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FastaIndexView):
            return False
        return list(self.values()) == list(other.values())

    def __repr__(self):
        return f"Fasta Index from {self._filename}\n" f"Full Header: {self._full_header}\n" f"Seqname: {self.keys()}"

    def __str__(self):
        return f"Fasta Index from {self._filename} [Full Header: {self._full_header}]"


def create_fai_from_fasta(fasta_filename: str, fai_filename: str, full_header: bool = True, show_tqdm: bool = True):
    """
    TODO: docs

     .. versionadded:: 1.0.2
    """
    FastaIndexView.from_fasta(filename=fasta_filename, full_header=full_header, show_tqdm=show_tqdm).write(fai_filename)
