"""
encode_fastq_to_parquet.py -- A TGS- and Big Data-Ready FASTQ to Apache Parquet Encoder
"""

import os
from typing import List

import numpy as np
import numpy.typing as npt
import pandas as pd

from labw_utils.bioutils.algorithm.sequence import decode_phred33
from labw_utils.bioutils.parser.fastq import FastqIterator
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def base_encoder(s: str) -> int:
    try:
        return "AGCT".index(s)
    except ValueError:
        return -1


def seq_encoder(seq: str, targeted_size: int) -> npt.NDArray:
    retn = np.zeros(targeted_size, dtype=np.float64)
    retn.fill(np.nan)
    retn[0:len(seq)] = np.array([
        base_encoder(s) for s in seq.upper()
    ], dtype=int)
    return retn


def qual_encoder(qual: str, targeted_size: int) -> npt.NDArray:
    retn = np.zeros(targeted_size, dtype=np.float64)
    retn.fill(np.nan)
    retn[0:len(qual)] = np.array(list(
        decode_phred33(qual)
    ), dtype=int)
    return retn


def encode(filepath: str):
    _lh.info("Start parsing '%s'...", filepath)
    if not os.path.exists(filepath):
        _lh.error("File '%s' not exist!", filepath)
    outdir_path = filepath + ".stats.d"
    os.makedirs(filepath + ".stats.d", exist_ok=True)
    with open(os.path.join(outdir_path, "all.tsv"), "wt") as all_writer:
        all_writer.write("\t".join((
            "SEQID",
            "GC",
            "LEN",
            "MEANQUAL"
        )) + "\n")
        max_length = 0
        num_records = 0
        for fastq_record in FastqIterator(filename=filepath):
            len_record = len(fastq_record)
            max_length = max(max_length, len_record)
            num_records += 1
        chunk_size = 10000
        col_names = list(str(n) for n in range(max_length))
        _lh.info("Would save encoded FASTQ information to %d * %d matrix", num_records, max_length)
        all_nt = None
        all_qual = None
        for i, fastq_record in enumerate(FastqIterator(filename=filepath)):
            j = i % chunk_size
            if j == 0:
                if all_nt is not None:
                    pd.DataFrame(
                        all_nt,
                        columns=col_names
                    ).to_parquet(f"{filepath}_nt.{i // chunk_size}.parquet")
                    pd.DataFrame(
                        all_qual,
                        columns=col_names
                    ).to_parquet(f"{filepath}_qual.{i // chunk_size}.parquet")
                all_nt = np.ndarray((chunk_size, max_length), dtype=np.float64)
                all_qual = np.ndarray((chunk_size, max_length), dtype=np.float64)
            all_nt[j, :] = seq_encoder(fastq_record.sequence, max_length)
            all_qual[j, :] = qual_encoder(fastq_record.quality, max_length)
        pd.DataFrame(
            all_nt,
            columns=col_names
        ).to_parquet(f"{filepath}_nt.{i // chunk_size}.parquet")
        pd.DataFrame(
            all_qual,
            columns=col_names
        ).to_parquet(f"{filepath}_qual.{i // chunk_size}.parquet")


def main(args: List[str]):
    for name in args:
        encode(name)
