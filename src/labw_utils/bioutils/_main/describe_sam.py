import json
import os
from typing import List, Optional, Union

from labw_utils import UnmetDependenciesError

try:
    import pytest

    pysam = pytest.importorskip("pysam")
except ImportError:
    try:
        import pysam
    except ImportError:
        raise UnmetDependenciesError("pysam")

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


def qc(
        sam_path: str
):
    if not os.path.exists(sam_path):
        _lh.error("Sam file at %s not exist!", sam_path)
        exit(1)
    out_dir_path = sam_path + ".stats.d"
    os.makedirs(out_dir_path, exist_ok=True)

    modestr = get_mode_str(sam_path=sam_path)

    with pysam.AlignmentFile(sam_path, modestr) as samfile:
        with open(os.path.join(out_dir_path, "header.json"), "w") as header_writer:
            header = samfile.header
            json.dump(header.to_dict(), header_writer)
    file_length = get_file_length(sam_path, modestr)
    determine_read_quality(
        sam_path=sam_path,
        modestr=modestr,
        out_dir_path=out_dir_path,
        file_length=file_length
    )


def main(args: List[str]):
    for name in args:
        qc(name)
