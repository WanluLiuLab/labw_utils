"""
Convert GTF to MTR -- Minimal Transcript Representation

We assume that the GTF was sorted in a form that exons are under corresponding transcript.

Two RMSK-TEs are removed at this step.
"""
import glob
import json
import multiprocessing
import os.path
import re
import uuid
import hashlib
from typing import Optional

import pandas as pd

from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.feature import FeatureType

UNASSIGNED_REGEX = re.compile(r"^unassigned_transcript_\d+$")


def convert_unassigned_transcript_id(tid: Optional[str]) -> str:
    if tid is None:
        return "unnamed_transcript_" + str(uuid.uuid4())
    if UNASSIGNED_REGEX.match(tid) is not None:
        return "unassigned_transcript_" + str(uuid.uuid4())
    return tid


def convert(
        src_gtf_file_path: str,
        dst_mtr_parquet_path: str
):
    mtr = {
        "SEQNAME": [],
        "START": [],
        "END": [],
        "STRAND": [],
        "TRANSCRIPT_ID_ORI": [],
        "TRANSCRIPT_ID": [],
        "FPATH": [],
        "FTSHA": [],
        "EXONS": [],
        "ATTRS": []
    }
    with GtfIterator(src_gtf_file_path, show_tqdm=True) as reader:
        for record in reader:
            if record.parsed_feature == FeatureType.TRANSCRIPT:
                mtr["SEQNAME"].append(record.seqname)
                mtr["START"].append(record.start)
                mtr["END"].append(record.end)
                mtr["STRAND"].append(record.strand)
                mtr["TRANSCRIPT_ID_ORI"].append(record.attribute_get("transcript_id", ""))
                new_transcript_id = convert_unassigned_transcript_id(record.attribute_get("transcript_id", None))
                mtr["TRANSCRIPT_ID"].append(new_transcript_id)
                mtr["FPATH"].append(src_gtf_file_path)
                mtr["FTSHA"].append(
                    hashlib.sha256(bytes(
                        src_gtf_file_path + new_transcript_id,
                        encoding="UTF-8"
                    )).hexdigest()
                )
                mtr["EXONS"].append([])
                mtr["ATTRS"].append(
                    json.dumps({k: v for k, v in zip(record.attribute_keys, record.attribute_values)})
                )
            else:
                assert record.seqname == mtr["SEQNAME"][-1]
                assert record.strand == mtr["STRAND"][-1]
                assert record.attribute_get("transcript_id", "") == mtr["TRANSCRIPT_ID_ORI"][-1]
                mtr["EXONS"][-1].append((record.start, record.end))
    mtr_df = pd.DataFrame(mtr)
    mtr_df.to_parquet(dst_mtr_parquet_path, index=False)


if __name__ == "__main__":
    os.makedirs("pre_processed_mtr", exist_ok=True)
    ppool = []
    for fn in glob.glob(os.path.join("pre_processed_gtf", "*.gtf")):
        dst_fn = os.path.join("pre_processed_mtr", os.path.basename(fn) + ".mtr.parquet")
        ppool.append(multiprocessing.Process(target=convert, args=(fn, dst_fn)))
        ppool[-1].start()
    for p in ppool:
        p.join()

    dfs = []
    for fn in glob.glob(os.path.join("pre_processed_mtr", "*.parquet")):
        dfs.append(pd.read_parquet(fn))
    pd.concat(dfs).to_parquet("final_mtr.parquet", index=False)
