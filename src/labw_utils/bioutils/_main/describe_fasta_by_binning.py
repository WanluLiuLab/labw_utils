"""
describe_fasta_by_binning.py -- Describe statistics on FASTA by binning.

.. versionadded:: 1.0.3
"""

from __future__ import annotations

__all__ = (
    "create_parser",
    "main",
)

import argparse
import json
import os.path
from collections import defaultdict

from labw_utils import UnmetDependenciesError
from labw_utils.bioutils.accession_matcher import infer_accession_type
from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Any, List, Dict

try:
    import numpy as np
    import numpy.typing as npt
except ImportError:
    raise UnmetDependenciesError("numpy")

try:
    import pandas as pd
except ImportError:
    raise UnmetDependenciesError("pandas")

try:
    import pyarrow

    DEST_FORMAT = "PARQUET"
except ImportError:
    DEST_FORMAT = "CSV"

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils describe_fasta_by_binning",
        description=__doc__.splitlines()[1],
    )
    parser = FrontendOptSpecs.patch(parser, "-f")
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output file basename. Will be {out}.csv or {out}.parquet if Apache Arrow is installed",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--nbins",
        required=False,
        help="Number of bins for each chromosome",
        nargs="?",
        type=int,
        action="store",
        default=400,
    )
    parser.add_argument(
        "--metadata_only",
        help="Dump metadata only",
        action="store_true",
    )
    return parser


def main(args: List[str]) -> None:
    out_metadata = {}
    args = create_parser().parse_args(args)
    nbins = args.nbins
    fasta_file_path = os.path.abspath(args.fasta)
    fa = FastaViewFactory(fasta_file_path, full_header=False, read_into_memory=False)
    out_metadata.update({"FATSA_FILE_PATH": fasta_file_path, "FASTA_CHRS": []})

    fail_to_idenfy = 0
    num_of_contigs = 0
    for chr_name in fa.chr_names:
        num_of_contigs += 1
        inf_type = infer_accession_type(chr_name)
        if inf_type is not None:
            inf_type = infer_accession_type(chr_name).as_dict()
        else:
            fail_to_idenfy += 1
        out_metadata["FASTA_CHRS"].append({"NAME": chr_name, "LEN": fa.get_chr_length(chr_name), "TYPE": inf_type})
    _lh.info(
        "Identified %d out of %d contigs (%.2f%%)",
        num_of_contigs - fail_to_idenfy,
        num_of_contigs,
        (num_of_contigs - fail_to_idenfy) / num_of_contigs * 100,
    )

    with get_writer(f"{args.out}.json") as metadata_writer:
        json.dump(out_metadata, metadata_writer, indent=4)
    if args.metadata_only:
        return
    out_dataframe: List[Dict[str, Any]] = []
    tqdm_total = sum(map(fa.get_chr_length, fa.chr_names))

    with tqdm(desc="Parsing FASTA", total=tqdm_total) as pbar:
        for chr_name in fa.chr_names:
            seq_len = fa.get_chr_length(chr_name)
            if seq_len < nbins:
                _lh.warning("Contig '%s' omitted: too short", chr_name)
            if seq_len < 20 * nbins:
                nbins = seq_len // 20
            segment_length = seq_len // nbins
            for start in range(0, seq_len - segment_length, segment_length):
                end = start + segment_length
                seq = fa.sequence(chr_name, start, end)
                nt_counts = defaultdict(lambda: 0)
                stats_dict = {"chr_name": chr_name, "start": start}
                for nt in seq:
                    nt_counts[nt] += 1
                stats_dict.update(nt_counts)
                out_dataframe.append(stats_dict)
                pbar.update(segment_length)
    _lh.info("Finished parsing, writing to disk...")
    out_dataframe_df = pd.DataFrame(out_dataframe)
    if DEST_FORMAT == "CSV":
        out_dataframe_df.to_csv(f"{args.out}.csv.xz", index=False)
    else:
        out_dataframe_df.to_parquet(f"{args.out}.parquet", index=False)
    _lh.info("Finished")
