import os
from collections import defaultdict
from typing import List, Any, Dict, Iterable

import pandas as pd

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.quantification_optimized_gene_tree import QuantificationOptimizedGeneTree
from labw_utils.bioutils.parser.feature import GtfIterator
from labw_utils.commonutils.importer.tqdm_importer import tqdm


def exclude_non_standard_chromosome_name(
        chr_names: Iterable[str]
) -> Iterable[str]:
    for chr_name in chr_names:
        if chr_name.startswith("N"):
            if not chr_name.startswith("NC_"):
                continue
        elif chr_name.startswith("chr"):
            if (
                    chr_name.startswith("chrUn") or
                    chr_name.endswith("_random") or
                    chr_name.endswith("_alt") or
                    chr_name.endswith("_decoy")
            ):
                continue
        else:
            if (
                    chr_name.startswith("KI") or
                    chr_name.startswith("GL") or
                    chr_name.startswith("KQ") or
                    chr_name.startswith("CHR_")
            ):
                continue
        yield chr_name


def main(args: List[str]) -> None:
    segment_length = 100000
    fasta_filename, gtf_filename = args
    fa = FastaViewFactory(fasta_filename, full_header=False, read_into_memory=True)
    gtf_intervals = QuantificationOptimizedGeneTree.from_feature_iterator(
        GtfIterator(gtf_filename),
        feature_attribute_name="transcript_id",
        feature_type="transcript"
    )
    summary_dirname = fasta_filename + ".summary.d"
    os.makedirs(summary_dirname, exist_ok=True)
    grf_chromosome_names = list(gtf_intervals.iter_chromosome_names())
    out_dataframe: List[Dict[str, Any]] = []
    for chr_name in exclude_non_standard_chromosome_name(fa.chr_names):
        if chr_name not in grf_chromosome_names:
            continue

        seq_len = fa.get_chr_length(chr_name)
        for start in tqdm(
                range(0, seq_len - segment_length, segment_length),
                desc=f"Iterating through {chr_name}"
        ):
            end = start + segment_length
            nt_counts = defaultdict(lambda: 0)
            nt_counts.update({i: 0 for i in "AGCTNagctn"})
            for nt in fa.sequence(chr_name, start, end):
                nt_counts[nt] += 1
            n_overlapping_gtf_intervals = sum(
                len(list(gtf_intervals.overlap(((chr_name, strand), start, end))))
                for strand in (True, False, None)
            )
            nt_counts.update({
                "chr_name": chr_name,
                "start": start,
                "gtf_intervals": n_overlapping_gtf_intervals
            })
            out_dataframe.append(nt_counts)
    out_dataframe_df = pd.DataFrame(out_dataframe)
    out_dataframe_df.to_csv(os.path.join(summary_dirname, "summary.csv.xz"), index=False)
