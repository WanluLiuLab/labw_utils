from typing import List

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory, split_fasta


def main(args: List[str]):
    for arg in args:
        split_fasta(FastaViewFactory(arg))
