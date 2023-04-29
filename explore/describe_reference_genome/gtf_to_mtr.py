"""
Convert GTF to MTR -- Minimal Transcript Representation

We assume that the GTF was sorted in a form that exons are under corresponding transcript.

Two RMSK-TEs are removed at this step.
"""
import glob
import hashlib
import json
import os.path
import re
import uuid

import pandas as pd

from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.record.feature import FeatureType
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Optional

UNASSIGNED_REGEX = re.compile(r"^unassigned_transcript_\d+$")
_lh = get_logger(__name__)
setup_basic_logger()


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
    fav = FastaViewFactory(
        os.path.join("pre_processed_fa", "hg38.fa"),
        read_into_memory=False
    )
    mtr = {
        "SEQNAME": [],
        "START": [],
        "END": [],
        "STRAND": [],
        "TRANSCRIPT_ID_ORI": [],
        "TRANSCRIPT_ID": [],
        "FPATH": [],
        "FTSHA": [],
        "EXON_BOUNDARIES": [],
        "SPLICE_SITES": [],
        "ATTRS": [],
        "CDNA": [],
        "CDNA_UNSPLICED": []
    }
    gt = DiploidGeneTree.from_gtf_file(
        src_gtf_file_path,
        show_tqdm=True,
        gene_implementation=DumbGene
    )
    for record in tqdm(gt.transcript_values, desc="Iteratong transcripts..."):
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
            mtr["EXON_BOUNDARIES"].append(
                list(record.exon_boundaries)
            )
            mtr["SPLICE_SITES"].append(
                list(record.splice_sites)
            )
            mtr["ATTRS"].append(
                json.dumps({k: v for k, v in zip(record.attribute_keys, record.attribute_values)})
            )
            mtr["CDNA"].append(
                record.transcribe(fav.sequence)
            )
            mtr["CDNA_UNSPLICED"].append(
                record.transcribe_unspliced(fav.sequence)
            )

    mtr_df = pd.DataFrame(mtr)
    mtr_df.to_parquet(dst_mtr_parquet_path, index=False)


if __name__ == "__main__":
    os.makedirs("pre_processed_mtr", exist_ok=True)
    for fn in glob.glob(os.path.join("pre_processed_gtf", "*.gtf")):
        dst_fn = os.path.join("pre_processed_mtr", os.path.basename(fn) + ".mtr.parquet")
        convert(fn, dst_fn)
