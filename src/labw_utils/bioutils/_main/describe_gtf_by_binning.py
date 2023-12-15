"""
describe_gtf_by_binning.py -- Describe number of features in GTF by binning.

.. versionadded:: 1.0.2

FIXME: This app uses GeneView v0.1!
"""
__all__ = (
    "create_parser",
    "main",
)

import argparse
import json
import os.path

from labw_utils import UnmetDependenciesError
from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.bioutils.datastructure.quantification_optimized_gene_tree import (
    QuantificationOptimizedGeneTree,
)
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper import pickle_helper
from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Any, Dict

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
_lh.warning("This module is not finished")  # TODO


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m labw_utils.bioutils describe_gtf_by_binning",
        description=__doc__.splitlines()[1],
    )
    parser = FrontendOptSpecs.patch(parser, "-g")
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
    return parser


def main(args: List[str]) -> None:
    out_metadata = {}
    args = create_parser().parse_args(args)
    fasta_file_path = os.path.abspath(args.fasta)
    fa = FastaViewFactory(fasta_file_path, full_header=False, read_into_memory=False)
    gtf_file_path = os.path.abspath(args.gtf)
    gv = GeneViewFactory.from_file(gtf_file_path)
    isoform_qoidx_path = gtf_file_path + ".isoform.qoidx.xz"
    gene_qoidx_path = gtf_file_path + ".gene.qoidx.xz"
    if os.path.exists(isoform_qoidx_path):
        gtf_isoform_intervals = pickle_helper.load(isoform_qoidx_path)
    else:
        gtf_isoform_intervals = QuantificationOptimizedGeneTree.from_feature_iterator(
            tqdm([gv.iter_genes()], desc="Creating Isoform-Level QOIDX"),
            feature_attribute_name="transcript_id",
            feature_type="transcript",
        )
        pickle_helper.dump(gtf_isoform_intervals, isoform_qoidx_path)
    if os.path.exists(gene_qoidx_path):
        gtf_gene_intervals = pickle_helper.load(gene_qoidx_path)
    else:
        gtf_gene_intervals = QuantificationOptimizedGeneTree.from_feature_iterator(
            tqdm([gv.iter_transcripts()], desc="Creating Gene-Level QOIDX"),
            feature_attribute_name="gene_id",
            feature_type="gene",
        )
        pickle_helper.dump(gtf_gene_intervals, gene_qoidx_path)
    fasta_chromosome_names = list(fa.chr_names)
    gtf_chromosome_names = gtf_isoform_intervals.iter_chromosome_names()
    final_chromosome_names = list(
        filter(
            lambda chr_name: chr_name in fasta_chromosome_names,
            gtf_isoform_intervals.iter_chromosome_names(),
        )
    )

    out_metadata.update(
        {
            "FATSA_FILE_PATH": fasta_file_path,
            "FASTA_CHR_NAME": fasta_chromosome_names,
            "FINAL_CHR_NAME": final_chromosome_names,
            "GTF_FILE_PATH": gtf_file_path,
            "GTF_CHR_NAME": gtf_chromosome_names,
        }
    )
    with get_writer(f"{args.out}.json") as metadata_writer:
        json.dump(out_metadata, metadata_writer)

    out_dataframe: List[Dict[str, Any]] = []
    tqdm_total = sum(map(fa.get_chr_length, final_chromosome_names))

    with tqdm(desc="Parsing GTF", total=tqdm_total) as pbar:
        for chr_name in final_chromosome_names:
            seq_len = fa.get_chr_length(chr_name)
            segment_length = seq_len // args.nbins
            for start in range(0, seq_len - segment_length, segment_length):
                end = start + segment_length
                stats_dict = {"chr_name": chr_name, "start": start}
                stats_dict.update(
                    {
                        "gtf_isoform_intervals_pos": len(
                            list(gtf_isoform_intervals.overlap(((chr_name, True), start, end)))
                        ),
                        "gtf_isoform_intervals_neg": len(
                            list(gtf_isoform_intervals.overlap(((chr_name, False), start, end)))
                        ),
                        "gtf_isoform_intervals_strandless": len(
                            list(gtf_isoform_intervals.overlap(((chr_name, None), start, end)))
                        ),
                        "gtf_gene_intervals_pos": len(list(gtf_gene_intervals.overlap(((chr_name, True), start, end)))),
                        "gtf_gene_intervals_neg": len(
                            list(gtf_gene_intervals.overlap(((chr_name, False), start, end)))
                        ),
                        "gtf_gene_intervals_strandless": len(
                            list(gtf_gene_intervals.overlap(((chr_name, None), start, end)))
                        ),
                    }
                )
                out_dataframe.append(stats_dict)
                pbar.update(segment_length)
    _lh.info("Finished parsing, writing to disk...")
    out_dataframe_df = pd.DataFrame(out_dataframe)
    if DEST_FORMAT == "CSV":
        out_dataframe_df.to_csv(f"{args.out}.csv.xz", index=False)
    else:
        out_dataframe_df.to_parquet(f"{args.out}.parquet", index=False)
    _lh.info("Finished")
