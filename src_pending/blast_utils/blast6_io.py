from __future__ import annotations

import pandas as pd

from blast_utils import NON_EXIST
from blast_utils.blast_db import BlastSearchConfig
from labw_utils.typing_importer import Mapping


def parse_blast6(
    src_blast6_df: pd.DataFrame,
    title_dict: Mapping[str, str],
    txid_dict: Mapping[str, int],
    seq_dict: Mapping[str, str],
    bsc: BlastSearchConfig,
):
    if len(src_blast6_df) == 0:
        return src_blast6_df
    src_blast6_df = (
        src_blast6_df.query(f"evalue < {bsc.filter_e_value}")
        .query(f"pident > {bsc.pident_cutoff}")
        .sort_values(by=["qseqid", "sseqid"])
    )

    if bsc.keep_longest_hit_only:
        src_blast6_df = (
            src_blast6_df.groupby(by=["qseqid", "sseqid"])
            .apply(lambda df: df.sort_values("length", ascending=False).head(1))
            .droplevel(0)
        )

    src_blast6_df["qstart"] = src_blast6_df["qstart"] - 1
    src_blast6_df["sstart"] = src_blast6_df["sstart"] - 1
    src_blast6_df["title"] = list(
        map(
            lambda x: title_dict.get(x, "UNKNOWN_TITLE"),
            src_blast6_df["sseqid"],
        )
    )
    src_blast6_df["txid"] = list(map(lambda x: txid_dict.get(x, NON_EXIST), src_blast6_df["sseqid"]))
    src_blast6_df["seq"] = list(
        map(
            lambda tup: seq_dict.get(f"{tup.sseqid}:{tup.sstart}-{tup.send}", "NNNN"),
            src_blast6_df.itertuples(index=False),
        )
    )
    return src_blast6_df


def read_blast6(blast_result_file: str) -> pd.DataFrame:
    return pd.read_csv(
        blast_result_file,
        sep="\t",
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
        engine="c",
        comment="#",
    )
