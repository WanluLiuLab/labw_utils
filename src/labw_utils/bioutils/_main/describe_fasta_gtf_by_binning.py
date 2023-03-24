import argparse
import json
import os.path
from collections import defaultdict
from typing import List, Any, Dict

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper import pickle_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

try:
    import pandas as pd
except ImportError:
    raise UnmetDependenciesError("pandas")

try:
    import pyarrow

    DEST_FORMAT = "PARQUET"
except ImportError:
    DEST_FORMAT = "CSV"

from labw_utils.bioutils.algorithm.sequence import is_valid_chrname
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.quantification_optimized_gene_tree import QuantificationOptimizedGeneTree
from labw_utils.bioutils.parser.feature import GtfIterator
from labw_utils.commonutils.importer.tqdm_importer import tqdm

_lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g", "--gtf",
        required=False,
        help="Gtf to filter from",
        nargs='?',
        type=str,
        action='store',
        default=None
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="Gtf to filter from",
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
    return parser.parse_args(args)


def main(args: List[str]) -> None:
    out_metadata = {}
    args = _parse_args(args)
    fasta_file_path = os.path.abspath(args.fasta)
    fa = FastaViewFactory(fasta_file_path, full_header=False, read_into_memory=False)
    if args.gtf is None:
        gtf_intervals = None
        gtf_file_path = None
        grf_chromosome_names = []
    else:
        gtf_file_path = os.path.abspath(args.gtf)
        qoidx_path = gtf_file_path + ".qoidx.xz"
        if os.path.exists(gtf_file_path + ".qoidx.xz"):
            gtf_intervals = pickle_helper.load(qoidx_path)
        else:
            gtf_intervals = QuantificationOptimizedGeneTree.from_feature_iterator(
                GtfIterator(gtf_file_path),
                feature_attribute_name="transcript_id",
                feature_type="transcript"
            )
            pickle_helper.dump(gtf_intervals, qoidx_path)
        grf_chromosome_names = list(gtf_intervals.iter_chromosome_names())

    fasta_chr_names = list(fa.chr_names)
    fasta_valid_chr_names = list(filter(is_valid_chrname, fasta_chr_names))
    out_metadata.update({
        "FATSA_FILE_PATH": fasta_file_path,
        "FASTA_CHR_NAME": fasta_valid_chr_names,
        "FASTA_VALID_CHR_NAME": fasta_valid_chr_names,
        "GTF_FILE_PATH": gtf_file_path,
        "GTF_CHR_NAME": grf_chromosome_names
    })

    if gtf_intervals is not None:
        fasta_valid_chr_names = list(filter(lambda _chr_name: _chr_name in grf_chromosome_names, fasta_valid_chr_names))

    out_dataframe: List[Dict[str, Any]] = []
    tqdm_total = sum(map(fa.get_chr_length, fasta_valid_chr_names))

    with tqdm(desc="Parsing FASTA", total=tqdm_total) as pbar:
        for chr_name in fasta_valid_chr_names:
            seq_len = fa.get_chr_length(chr_name)
            segment_length = seq_len // args.nbins
            for start in range(0, seq_len - segment_length, segment_length):
                end = start + segment_length
                nt_counts = defaultdict(lambda: 0)
                nt_counts.update({
                    "chr_name": chr_name,
                    "start": start,
                })
                nt_counts.update({i: 0 for i in "AGCTNagctn"})
                for nt in fa.sequence(chr_name, start, end):
                    nt_counts[nt] += 1
                if gtf_intervals is not None:
                    gtf_intervals_pos = len(list(gtf_intervals.overlap(((chr_name, True), start, end))))
                    gtf_intervals_neg = len(list(gtf_intervals.overlap(((chr_name, False), start, end))))
                    gtf_intervals_strandless = len(list(gtf_intervals.overlap(((chr_name, None), start, end))))
                else:
                    gtf_intervals_pos = None
                    gtf_intervals_neg = None
                    gtf_intervals_strandless = None
                nt_counts.update({
                    "gtf_intervals_pos": gtf_intervals_pos,
                    "gtf_intervals_neg": gtf_intervals_neg,
                    "gtf_intervals_strandless": gtf_intervals_strandless,
                })
                out_dataframe.append(nt_counts)
                pbar.update(segment_length)
    _lh.info("Finished parsing, writing to disk...")
    out_dataframe_df = pd.DataFrame(out_dataframe)
    if DEST_FORMAT == "CSV":
        out_dataframe_df.to_csv(f"{args.out}.csv.xz", index=False)
    else:
        out_dataframe_df.to_parquet(f"{args.out}.parquet", index=False)
    with get_writer(f"{args.out}.json") as metadata_writer:
        json.dump(out_metadata, metadata_writer)
    _lh.info("Finished")
