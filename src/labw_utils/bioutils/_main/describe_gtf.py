"""
describe_gtf.py -- Get statistics about GTF files that can be parsed into a Gene-Transcript-Exon Three-Tier Structure

SYNOPSIS: python -m labw_utils.bioutils describe_gtf [GTF] [[GTF]...]

where [GTF] are path to GTF files you wish to describe.
"""

__all__ = (
    "main",
)

from typing import List

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_helper import describe


def main(args: List[str]):
    for arg in args:
        describe(arg, arg)
