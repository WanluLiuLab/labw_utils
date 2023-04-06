import argparse
import json
import os.path
from typing import List, Any, Dict

from labw_utils import UnmetDependenciesError
from labw_utils.bioutils.algorithm.sequence import is_valid_chrname
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.bioutils.datastructure.quantification_optimized_gene_tree import QuantificationOptimizedGeneTree
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper import pickle_helper
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
        "-g", "--gtf",
        required=True,
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


def sdi(abundance_data: npt.ArrayLike) -> float:
    return 1 - np.sum(np.power(abundance_data / np.sum(abundance_data), 2))


def main(args: List[str]) -> None:
    out_metadata = {}
    args = _parse_args(args)
    fasta_file_path = os.path.abspath(args.fasta)
    fa = FastaViewFactory(fasta_file_path, full_header=False, read_into_memory=False)
    gtf_file_path = os.path.abspath(args.gtf)
    gv = GeneViewFactory.from_file(gtf_file_path)
    isoform_qoidx_path = gtf_file_path + ".isoform.qoidx.xz"
    gene_qoidx_path = gtf_file_path + ".gene.qoidx.xz"
    if os.path.exists(isoform_qoidx_path):
        gtf_isoform_intervals = pickle_helper.load(isoform_qoidx_path)
    else:
        gtf_isoform_intervals = None
    if os.path.exists(gene_qoidx_path):
        gtf_gene_intervals = pickle_helper.load(gene_qoidx_path)
    else:
        gtf_gene_intervals = None
    if gtf_isoform_intervals is None or gtf_gene_intervals is None:
        if gtf_isoform_intervals is None:
            gtf_isoform_intervals = QuantificationOptimizedGeneTree.from_feature_iterator(
                tqdm([gv.iter_genes()], desc="Creating Isoform-Level QOIDX"),
                feature_attribute_name="transcript_id",
                feature_type="transcript"
            )
            pickle_helper.dump(gtf_isoform_intervals, isoform_qoidx_path)
        if gtf_gene_intervals is None:
            gtf_gene_intervals = QuantificationOptimizedGeneTree.from_feature_iterator(
                tqdm([gv.iter_transcripts()], desc="Creating Gene-Level QOIDX"),
                feature_attribute_name="gene_id",
                feature_type="gene"
            )
            pickle_helper.dump(gtf_gene_intervals, gene_qoidx_path)
    grf_chromosome_names = list(gtf_isoform_intervals.iter_chromosome_names())

    fasta_chr_names = list(fa.chr_names)
    fasta_chr_lengths = list(fa.get_chr_length(chr_name) for chr_name in fasta_chr_names)
    fasta_valid_chr_names = list(filter(is_valid_chrname, fasta_chr_names))
    out_metadata.update({
        "FATSA_FILE_PATH": fasta_file_path,
        "FASTA_CHR_NAME": fasta_valid_chr_names,
        "FASTA_CHR_LENGTH": fasta_chr_lengths,
        "FASTA_VALID_CHR_NAME": fasta_valid_chr_names,
        "GTF_FILE_PATH": gtf_file_path,
        "GTF_CHR_NAME": grf_chromosome_names
    })
    with get_writer(f"{args.out}.json") as metadata_writer:
        json.dump(out_metadata, metadata_writer)
    if gtf_isoform_intervals is not None:
        fasta_valid_chr_names = list(filter(lambda _chr_name: _chr_name in grf_chromosome_names, fasta_valid_chr_names))

    out_dataframe: List[Dict[str, Any]] = []
    tqdm_total = sum(map(fa.get_chr_length, fasta_valid_chr_names))

    with tqdm(desc="Parsing FASTA", total=tqdm_total) as pbar:
        for chr_name in grf_chromosome_names:
            try:
                seq_len = fa.get_chr_length(chr_name)
            except KeyError:
                continue  # pass
            segment_length = seq_len // args.nbins
            for start in range(0, seq_len - segment_length, segment_length):
                end = start + segment_length
                stats_dict = {
                    "chr_name": chr_name,
                    "start": start
                }
                gtf_isoform_intervals_pos = len(list(gtf_isoform_intervals.overlap(((chr_name, True), start, end))))
                gtf_isoform_intervals_neg = len(list(gtf_isoform_intervals.overlap(((chr_name, False), start, end))))
                gtf_isoform_intervals_strandless = len(
                    list(gtf_isoform_intervals.overlap(((chr_name, None), start, end))))
                gtf_gene_intervals_pos = len(list(gtf_gene_intervals.overlap(((chr_name, True), start, end))))
                gtf_gene_intervals_neg = len(list(gtf_gene_intervals.overlap(((chr_name, False), start, end))))
                gtf_gene_intervals_strandless = len(list(gtf_gene_intervals.overlap(((chr_name, None), start, end))))

                stats_dict.update({
                    "gtf_isoform_intervals_pos": gtf_isoform_intervals_pos,
                    "gtf_isoform_intervals_neg": gtf_isoform_intervals_neg,
                    "gtf_isoform_intervals_strandless": gtf_isoform_intervals_strandless,
                    "gtf_gene_intervals_pos": gtf_gene_intervals_pos,
                    "gtf_gene_intervals_neg": gtf_gene_intervals_neg,
                    "gtf_gene_intervals_strandless": gtf_gene_intervals_strandless
                })
                out_dataframe.append(stats_dict)
                pbar.update(segment_length)
    _lh.info("Finished parsing, writing to disk...")
    out_dataframe_df = pd.DataFrame(out_dataframe)
    if DEST_FORMAT == "CSV":
        out_dataframe_df.to_csv(f"{args.out}.csv.xz", index=False)
    else:
        out_dataframe_df.to_parquet(f"{args.out}.parquet", index=False)
    _lh.info("Finished")
