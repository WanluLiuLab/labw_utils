"""
split_fasta.py -- Split input FASTA file into one-line FASTAs with one file per contig.

SYNOPSIS: python -m labw_utils.bioutils split_fasta [FASTA] [[FASTA]...]

where [FASTA] are path to FASTA files.

.. versionadded:: 1.0.2
"""

__all__ = "main"

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory, split_fasta
from labw_utils.typing_importer import List


def main(args: List[str]):
    if "--help" in args or "-h" in args:
        print(__doc__)
        return 0
    for arg in args:
        split_fasta(FastaViewFactory(arg))
