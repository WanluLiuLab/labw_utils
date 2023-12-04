"""
normalize_gtf.py -- Performs GTF normalization, etc.

.. versionadded:: 1.0.2
"""

__all__ = (
    "create_parser",
    "main",
)

import argparse

from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_feature_proxy import (
    DEFAULT_SORT_EXON_EXON_STRAND_POLICY,
    VALID_SORT_EXON_EXON_STRAND_POLICY,
)
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.old_feature_parser import GtfIterator, GtfWriter
from labw_utils.bioutils.record.feature import VALID_GTF_QUOTE_OPTIONS, DEFAULT_GTF_QUOTE_OPTIONS
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils normalize_gtf",
        description=__doc__.splitlines()[1],
    )
    parser = FrontendOptSpecs.patch(parser, "-g")
    parser.add_argument(
        "--three_tier",
        help="Whether to parse the GTF into Gene-Transcript-Exon Three-Tier Structure."
        + "Other features will be discarded, and missing genes & transcripts will be added by maximum span length",
        action="store_true",
    )
    parser.add_argument(
        "--quote",
        required=False,
        help="Whether to add quotes in alternative field of output GTF",
        nargs="?",
        type=str,
        action="store",
        choices=VALID_GTF_QUOTE_OPTIONS,  # TODO: Change to Enum
        default=DEFAULT_GTF_QUOTE_OPTIONS,
    )
    parser.add_argument(
        "--sort_exon_exon_number_policy",
        required=False,
        help="How you would rank exon_number? On Reference: unstranded. Sorted by bedtools: stranded. Need --three_tier",
        nargs="?",
        type=str,
        action="store",
        choices=VALID_SORT_EXON_EXON_STRAND_POLICY,  # TODO: Change to Enum
        default=DEFAULT_SORT_EXON_EXON_STRAND_POLICY,
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path to normalized output",
        nargs="?",
        type=str,
        action="store",
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    if args.three_tier:
        gv = GeneViewFactory.from_file(args.gtf, not_save_index=True)
        gv.standardize(sort_exon_exon_number_policy=args.sort_exon_exon_number_policy)
        gi = gv.get_iterator()
    else:
        gi = GtfIterator(args.gtf)
    with GtfWriter(args.out) as gw:
        for feature in gi:
            gw.write_feature(feature, quote=args.quote)
