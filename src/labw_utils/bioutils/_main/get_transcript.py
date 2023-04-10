"""
get_transcript.py -- Filter GTF records by a specific attributes
"""

__all__ = (
    "create_parser",
    "main"
)

import argparse
from typing import List

from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_helper import subset_gtf_by_attribute_value
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m labw_utils.bioutils get_transcript",
        description=__doc__.splitlines()[1]
    )
    parser = FrontendOptSpecs.patch(parser, "-g")
    parser.add_argument(
        "--attribute_name",
        required=False,
        help="Attribute to be filtered",
        nargs='?',
        type=str,
        action='store',
        default="transcript_id"
    )
    parser.add_argument(
        "--attribute_values",
        required=True,
        help="File path to file that contains legal values, once per line, can be quoted",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Filtered output GTF",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    possible_values = []
    for line in get_tqdm_line_reader(args.field_value):
        line = line.strip().strip("\"\'")  # Get rid of quotation marks produced by R
        possible_values.append(line)
    lh.info(f"{len(possible_values)} values loaded")
    subset_gtf_by_attribute_value(
        attribute_values=iter(possible_values),
        attribute_name=args.field_name,
        gtf_filename=args.gtf,
        out_filename=args.out
    )
