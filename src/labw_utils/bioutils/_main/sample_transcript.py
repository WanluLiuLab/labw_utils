"""
sample_transcript.py -- Sample fraction of transcripts inside GTF file.

.. versionadded:: 1.0.2
"""

__all__ = (
    "create_parser",
    "main",
)

import argparse
import random

from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.datastructure.gene_tree_helper import subset_gtf_by_attribute_value
from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils sample_transcript",
        description=__doc__.splitlines()[1],
    )
    parser = FrontendOptSpecs.patch(parser, "-g")
    parser.add_argument(
        "--percent",
        required=False,
        help="How many percent of transcript to be sampled",
        nargs="?",
        type=float,
        action="store",
        default=50,
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path to filtered output",
        nargs="?",
        type=str,
        action="store",
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    gi = GtfIterator(args.gtf)
    transcript_ids = set()
    for gtf_record in gi:
        tmp_tid = gtf_record.attribute_get("transcript_id", None)
        if tmp_tid is not None:
            transcript_ids.add(tmp_tid)
    transcript_ids = random.sample(list(transcript_ids), int(len(transcript_ids) * args.percent / 100))

    subset_gtf_by_attribute_value(
        attribute_values=transcript_ids,
        attribute_name="transcript_id",
        gtf_filename=args.gtf,
        out_filename=args.out,
    )
