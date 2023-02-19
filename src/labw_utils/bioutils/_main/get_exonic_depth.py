"""
get_transcript.py -- Filter GTF records by a specific attributes
"""

import argparse
import itertools
import os
from collections import defaultdict
from typing import List

import pysam

from labw_utils.bioutils.datastructure.gv_feature_proxy import merge_intervals
from labw_utils.bioutils._main.describe_sam import get_mode_str, get_file_length, turn_none_to_zero
from labw_utils.bioutils.parser.feature import GtfIterator
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g", "--gtf",
        required=False,
        help="Gtf to filter from",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        "-s", "--sam_path",
        required=True,
        help="SAM/BAM file to get depth from",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        "-l", "--known_mappable_length",
        required=False,
        help="If the mappable length is known, put it here",
        nargs='?',
        type=int,
        action='store',
        default=0
    )
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    sam_path = args.sam_path
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
    print(f"Totally primarily mapped {primary_mapped_bases} bases")
    if args.known_mappable_length == 0:
        gtf_intervals = defaultdict(lambda: [])
        for gtf_record in GtfIterator(args.gtf):
            if gtf_record.feature == "exon":
                gtf_intervals[gtf_record.seqname].append([gtf_record.start, gtf_record.end - 1])
        gtf_intervals = {
            k: merge_intervals(gtf_intervals[k]) for k in gtf_intervals
        }
        gtf_mappable_length = sum(
            interval[1] - interval[0] for interval in itertools.chain(*gtf_intervals.values())
        )
        print(f"GTF mappable length: {gtf_mappable_length}")
    else:
        gtf_mappable_length = args.known_mappable_length
    print(f"Depth: {primary_mapped_bases / gtf_mappable_length}")
