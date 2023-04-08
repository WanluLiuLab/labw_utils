"""
get_transcript.py -- Filter GTF records by a specific attributes
"""

import argparse
import itertools
import os
from collections import defaultdict
from typing import List, Optional, Union

from labw_utils import UnmetDependenciesError

try:
    import pytest

    pysam = pytest.importorskip("pysam")
except ImportError:
    pytest = None
    try:
        import pysam
    except ImportError:
        raise UnmetDependenciesError("pysam")

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_feature_proxy import merge_intervals
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.old_feature_parser import GtfIterator
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def get_file_length(
        sam_path: str,
        modestr: str
) -> int:
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


def determine_read_quality(
        sam_path: str,
        modestr: str,
        out_dir_path: str,
        file_length: int
):
    _lh.info("Determining read quality...")
    with pysam.AlignmentFile(sam_path, modestr) as samfile:
        with open(os.path.join(out_dir_path, "read_stat.tsv"), "w") as read_stat_writer:
            read_stat_writer.write("\t".join((
                "QUERY_NAME",
                "MAP_STAT",
                "QUERY_LENGTH",
                "REFERENCE_LENGTH",
                "CIGAR_INFERRED_QUERY_LENGTH",
                "CIGAR_INFERRED_READ_LENGTH",
                "MAPPING_QUALITY"
            )) + "\n")
            for read in tqdm(samfile.fetch(), total=file_length):
                read: pysam.AlignedSegment
                if read.is_unmapped:
                    map_stat = "unmapped"
                elif read.is_supplementary:
                    map_stat = "supplementary"
                elif read.is_secondary:
                    map_stat = "secondary"
                else:
                    map_stat = "primiary"
                read_stat_writer.write("\t".join((
                    read.query_name,  # "QUERY_NAME",
                    map_stat,  # "MAP_STAT",
                    str(turn_none_to_zero(read.query_length)),  # "QUERY_LENGTH",
                    str(turn_none_to_zero(read.reference_length)),  # "REFERENCE_LENGTH",
                    str(turn_none_to_zero(read.infer_query_length())),  # "CIGAR_INFERRED_QUERY_LENGTH",
                    str(turn_none_to_zero(read.infer_read_length())),  # "CIGAR_INFERRED_READ_LENGTH",
                    str(turn_none_to_zero(read.mapping_quality))  # "MAPPING_QUALITY"
                )) + "\n")


def get_mode_str(sam_path: str) -> str:
    if sam_path.endswith(".sam"):
        modestr = "r"
    elif sam_path.endswith(".bam"):
        modestr = "rb"
    else:
        _lh.error("Sam file at %s have unknown extensions!", sam_path)
        exit(1)
    return modestr

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
