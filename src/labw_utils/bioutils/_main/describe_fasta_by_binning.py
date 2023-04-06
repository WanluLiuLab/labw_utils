import argparse
import json
import os.path
from collections import defaultdict
from typing import List, Any, Dict

from labw_utils import UnmetDependenciesError
from labw_utils.bioutils.accession_matcher.match import infer_accession_type
from labw_utils.bioutils.algorithm.sequence import is_valid_chrname, get_gc_percent
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

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


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="FASTA to filter from",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        "-o", "--out",
        required=True,
        help="Output file basename. Will be {out}.csv or {out}.parquet if Apache Arrow is installed",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        "--nbins",
        required=False,
        help="Number of bins for each chromosome",
        nargs='?',
        type=int,
        action='store',
        default=400
    )
    parser.add_argument(
        "--metadata_only",
        help="Dump metadata only",
        action='store_true'
    )
    return parser.parse_args(args)


def sdi(abundance_data: npt.ArrayLike) -> float:
    return 1 - np.sum(np.power(abundance_data / np.sum(abundance_data), 2))


def main(args: List[str]) -> None:
    out_metadata = {}
    args = _parse_args(args)
    nbins = args.nbins
    fasta_file_path = os.path.abspath(args.fasta)
    fa = FastaViewFactory(fasta_file_path, full_header=False, read_into_memory=False)
    out_metadata.update({
        "FATSA_FILE_PATH": fasta_file_path,
        "FASTA_CHRS":[]
    })

    for chr_name in fa.chr_names:
        inf_type = infer_accession_type(chr_name)
        if inf_type is not None:
            inf_type = infer_accession_type(chr_name).as_dict()
        out_metadata["FASTA_CHRS"].append({
            "NAME": chr_name,
            "LEN": fa.get_chr_length(chr_name),
            "TYPE": inf_type
        })

    with get_writer(f"{args.out}.json") as metadata_writer:
        json.dump(out_metadata, metadata_writer)
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
                nbins = seq_len / 20
            segment_length = seq_len // nbins
            for start in range(0, seq_len - segment_length, segment_length):
                end = start + segment_length
                seq = fa.sequence(chr_name, start, end)
                nt_counts = defaultdict(lambda: 0)
                stats_dict = {
                    "chr_name": chr_name,
                    "start": start,
                    "gc": get_gc_percent(seq)
                }
                for nt in seq:
                    nt_counts[nt] += 1
                stats_dict["sdi"] = sdi(list(nt_counts.values()))
                out_dataframe.append(stats_dict)
                pbar.update(segment_length)
    _lh.info("Finished parsing, writing to disk...")
    out_dataframe_df = pd.DataFrame(out_dataframe)
    if DEST_FORMAT == "CSV":
        out_dataframe_df.to_csv(f"{args.out}.csv.xz", index=False)
    else:
        out_dataframe_df.to_parquet(f"{args.out}.parquet", index=False)
    _lh.info("Finished")
