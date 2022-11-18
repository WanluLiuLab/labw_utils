import os
from collections import defaultdict
from typing import List

import matplotlib.pyplot as plt
import pandas as pd

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm


def main(args: List[str]) -> None:
    plot_kwargs = dict(kind='area', stacked=True, figsize=(20, 10))
    segment_length = 1000000
    for fasta_filename in args:
        fa = FastaViewFactory(fasta_filename, full_header=False, read_into_memory=False)
        summary_dirname = fasta_filename + ".summary.d"
        os.makedirs(summary_dirname, exist_ok=True)
        for chr_name in fa.chr_names:
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
                        chr_name.startswith("KQ")
                ):
                    continue

            seq_len = fa.get_chr_length(chr_name)
            nt_counts_across_bins = []
            cat_counts_across_bins = []
            this_segment_length = segment_length
            if seq_len < 100:
                continue
            else:
                while seq_len // this_segment_length < 100:
                    this_segment_length //= 10
            for i in tqdm(
                    range(0, seq_len - this_segment_length, this_segment_length),
                    desc=f"Iterating through {chr_name}"
            ):
                seq = fa.sequence(chr_name, i, i + this_segment_length)
                nt_counts = defaultdict(lambda: 0)
                nt_counts.update({i: 0 for i in "AGCTNagct"})
                for nt in seq:
                    nt_counts[nt] += 1
                nt_counts_across_bins.append(nt_counts)
                cat_counts = {
                    "Soft Masked": sum(nt_counts[i] for i in "agct"),
                    "Hard Masked": sum(nt_counts[i] for i in "N"),
                    "Normal": sum(nt_counts[i] for i in "AGCT"),
                }
                cat_counts_across_bins.append(cat_counts)
            full_nt_counts = pd.DataFrame(nt_counts_across_bins).fillna(0).astype(int)
            full_nt_counts.plot(**plot_kwargs)
            plt.savefig(
                os.path.join(summary_dirname, chr_name + "_c_nt.png")
            )
            plt.close()
            full_cat_counts = pd.DataFrame(cat_counts_across_bins).fillna(0).astype(int)
            full_cat_counts.plot(**plot_kwargs)
            plt.savefig(
                os.path.join(summary_dirname, chr_name + "_c_cat.png")
            )
            plt.close()
