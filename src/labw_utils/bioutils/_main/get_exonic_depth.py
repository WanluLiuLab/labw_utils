"""
get_exonic_depth.py -- Get depth from RNA-Seq alignment files.

.. versionadded:: 1.0.2
"""

__all__ = (
    "create_parser",
    "main",
)

import argparse
import itertools
import os
from collections import defaultdict

from labw_utils import UnmetDependenciesError
from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.feature import FeatureType
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List, Optional, Union, Literal

try:
    import pytest

    pysam = pytest.importorskip("pysam")
except ImportError:
    pytest = None
    try:
        import pysam
    except ImportError as e:
        raise UnmetDependenciesError("pysam") from e

from labw_utils.bioutils.algorithm.utils import merge_intervals
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils get_exonic_depth",
        description=__doc__.splitlines()[1],
    )
    parser = FrontendOptSpecs.patch(parser, "-s")
    meg = parser.add_mutually_exclusive_group()
    FrontendOptSpecs.patch(meg, "-g", required=False)
    meg.add_argument(
        "-l",
        "--known_mappable_length",
        required=False,
        help="If the mappable length is known, put it here",
        nargs="?",
        type=int,
        action="store",
        default=0,
    )
    return parser


def get_file_length(sam_path: str, modestr: Literal["r", "rb"]) -> int:
    _lh.info("Determining file length...")
    file_length = 0
    with pysam.AlignmentFile(sam_path, modestr) as samfile:
        for _ in samfile.fetch():
            file_length += 1
    return file_length


def turn_none_to_zero(i: Optional[Union[int, float]]) -> Union[int, float]:
    if i is None:
        i = 0
    return i


def get_mode_str(sam_path: str) -> str:
    if sam_path.endswith(".sam"):
        modestr = "r"
    elif sam_path.endswith(".bam"):
        modestr = "rb"
    else:
        _lh.error("Sam file at %s have unknown extensions!", sam_path)
        exit(1)
    return modestr


def main(args: List[str]):
    args = create_parser().parse_args(args)
    sam_path = args.sam
    if not os.path.exists(sam_path):
        _lh.error("Sam file at %s not exist!", sam_path)
        exit(1)
    modestr = get_mode_str(sam_path=sam_path)
    file_length = get_file_length(sam_path, modestr)
    primary_mapped_bases = 0
    with pysam.AlignmentFile(sam_path, modestr) as samfile:
        for read in tqdm(samfile.fetch(), total=file_length):
            read: pysam.AlignedSegment
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            primary_mapped_bases += turn_none_to_zero(read.infer_query_length())
    _lh.info(f"Totally primarily mapped %d bases", primary_mapped_bases)
    if args.known_mappable_length == 0:
        _lh.info(f"Iterating intervals from GTF...")
        gtf_intervals = defaultdict(lambda: [])
        for gtf_record in GtfIterator(args.gtf):
            if gtf_record.parsed_feature == FeatureType.EXON:
                gtf_intervals[gtf_record.seqname].append([gtf_record.start, gtf_record.end - 1])
        _lh.info(f"Merging intervals...")
        gtf_intervals = {k: merge_intervals(gtf_intervals[k]) for k in gtf_intervals}
        gtf_mappable_length = sum(interval[1] - interval[0] for interval in itertools.chain(*gtf_intervals.values()))
        _lh.info(f"GTF mappable length: %d", gtf_mappable_length)
    else:
        gtf_mappable_length = args.known_mappable_length
    _lh.info(f"Depth: %.2f", primary_mapped_bases / gtf_mappable_length)
