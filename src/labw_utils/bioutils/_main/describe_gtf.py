"""
describe_gtf.py -- Get statistics about GTF files that can be parsed into a Gene-Transcript-Exon Three-Tier Structure

SYNOPSIS: python -m labw_utils.bioutils describe_gtf [GTF] [[GTF]...]

where [GTF] are path to GTF files you wish to describe.

.. versionadded:: 1.0.2

FIXME: This app uses GeneView v0.1!
"""

__all__ = ("main",)

from labw_utils.bioutils.datastructure.gene_tree_helper import describe
from labw_utils.typing_importer import List


def main(args: List[str]):
    for arg in args:
        describe(arg, arg)
