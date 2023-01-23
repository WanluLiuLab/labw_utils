import json
import os
import shutil
import subprocess
import threading
from typing import List

import pysam
import tqdm

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
            for read in tqdm.tqdm(samfile.fetch(), total=file_length):
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
                    str(read.query_length),  # "QUERY_LENGTH",
                    str(read.reference_length),  # "REFERENCE_LENGTH",
                    str(read.infer_query_length()),  # "CIGAR_INFERRED_QUERY_LENGTH",
                    str(read.infer_read_length()),  # "CIGAR_INFERRED_READ_LENGTH",
                    str(read.mapping_quality)  # "MAPPING_QUALITY"
                )) + "\n")


def get_mode_str(sam_path: str) -> str:
    if sam_path.endswith(".sam"):
        modestr = "r"
    elif sam_path.endswith(".bam"):
        modestr = "rb"
    else:
        _lh.error(f"Sam file at %s have unknown extensions!", sam_path)
        exit(1)
    return modestr


def determine_pileup_quality_lightweighted(
        sam_path: str,
        modestr: str,
        out_dir_path: str,
        num_all_pos: int
):
    _ = modestr, num_all_pos
    del modestr, num_all_pos
    _lh.info("Determining pileup quality...")
    with open(os.path.join(out_dir_path, "pileup_stat.tsv.gz"), "wb") as pileup_stat_writer:
        gz_p = subprocess.Popen(
            [
                "pigz",
                "-9cf",
                "-"
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE
        )
        gz_p.stdin.write(bytes("\t".join((
            "REFERENCE_NAME",
            "REFERENCE_POS",
            "NUM_READS"
        )) + "\n", encoding="utf-8"))
        samtools_p = subprocess.Popen(
            ["samtools", "depth", sam_path],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE
        )
        t2 = threading.Thread(target=shutil.copyfileobj, args=(gz_p.stdout, pileup_stat_writer))
        t2.start()
        t1 = threading.Thread(target=shutil.copyfileobj, args=(samtools_p.stdout, gz_p.stdin))
        t1.start()
        t1.join()
        gz_p.stdin.close()
        t2.join()
        samtools_p.wait()
        gz_p.wait()


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
            num_all_pos = sum(header.lengths)

    determine_pileup_quality_lightweighted(
        sam_path=sam_path,
        modestr=modestr,
        out_dir_path=out_dir_path,
        num_all_pos=num_all_pos
    )
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
