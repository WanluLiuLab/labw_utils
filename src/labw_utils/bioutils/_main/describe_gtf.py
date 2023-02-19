"""Get statistics about GTF files that can be parsed into a Gene-Transcript-Exon Three-Tier Structure"""

import statistics
from typing import List

from labw_utils.bioutils.datastructure.gene_view import GeneViewFactory
<<<<<<<< HEAD:src/labw_utils/bioutils/_main/get_gtf_statistics.py
from labw_utils.bioutils.parser.feature import GtfWriter
from matplotlib import pyplot as plt

========
>>>>>>>> origin/0.1.X-dev:src/labw_utils/bioutils/_main/describe_gtf.py
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer


def stat(item: List[int], fig_name: str):
    plt.clf()
    quantiles = statistics.quantiles(item)
    print(f"{fig_name} " + ", ".join((
        f"min={min(item)}",
        f"mean={round(statistics.mean(item), 2)}",
        f"median={statistics.median(item)}",
        f"max={max(item)}",
        f"quantiles={quantiles}"
    )))
    plt.hist(item, bins=150)
    # plt.xlim(0, quantiles[2])
    plt.savefig(f"{fig_name}_distribution.png")


def describe(input_filename: str, out_basename: str):
    gv = GeneViewFactory.from_file(input_filename, not_save_index=True)

    with get_writer(f"{out_basename}.gene.tsv") as gene_writer, \
            get_writer(f"{out_basename}.transcripts.tsv") as transcripts_writer, \
            get_writer(f"{out_basename}.exons.tsv") as exons_writer:
        gene_writer.write("\t".join((
            "GENE_ID",
            "TRANSCRIPT_NUMBER",
            "NAIVE_LENGTH",
            "TRANSCRIBED_LENGTH",
            "MAPPABLE_LENGTH"
        )) + "\n")
        transcripts_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "GENE_ID",
            "NAIVE_LENGTH",
            "TRANSCRIBED_LENGTH",
            "EXON_NUMBER"
        )) + "\n")
        exons_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "EXON_NUMBER",
            "NAIVE_LENGTH"
        )) + "\n")

        for gene in tqdm(desc="Iterating over genes...", iterable=gv.iter_genes()):

            gene_writer.write("\t".join((
                str(gene.gene_id),
                str(gene.number_of_transcripts),
                str(gene.naive_length),
                str(gene.transcribed_length),
                str(gene.mappable_length)
            )) + "\n")

            transcripts = list(gene.iter_transcripts())
            for t_i in range(len(transcripts)):
                transcript = transcripts[t_i]
                for exon in list(transcript.iter_exons()):
                    exons_writer.write("\t".join((
                        exon.transcript_id,
                        str(exon.exon_number),
                        str(exon.naive_length)
                    )) + "\n")
                transcripts_writer.write("\t".join((
                    transcript.transcript_id,
                    transcript.gene_id,
                    str(transcript.naive_length),
                    str(transcript.transcribed_length),
                    str(transcript.number_of_exons)
                )) + "\n")


def main(args: List[str]):
    for arg in args:
        describe(arg, arg)
