"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

import os
import re

from labw_utils.bioutils.algorithm.sequence import get_gc_percent
from labw_utils.bioutils.datastructure.fasta_view import FastaViewType
from labw_utils.bioutils.datastructure.gene_tree import GeneTreeInterface, GeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene, Gene
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.parser.gtf import GtfIterator, GtfIteratorWriter
from labw_utils.commonutils.appender import load_table_appender_class, TableAppenderConfig
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Iterable, Sequence, List

_lh = get_logger(__name__)


def transcribe_transcripts(
    it: Iterable[Transcript],
    dst_fasta_path: str,
    fv: FastaViewType,
    write_single_transcript: bool = True,
):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    intermediate_fasta_dir = ""
    if write_single_transcript:
        intermediate_fasta_dir = dst_fasta_path + ".d"
        os.makedirs(intermediate_fasta_dir, exist_ok=True)
    with get_writer(dst_fasta_path) as fasta_writer, load_table_appender_class("TSVTableAppender")(
        dst_fasta_path + ".stats",
        (
            "TRANSCRIPT_ID",
            "GENE_ID",
            "SEQNAME",
            "START",
            "END",
            "STRAND",
            "ABSOLUTE_LENGTH",
            "TRANSCRIBED_LENGTH",
            "GC",
        ),
        tac=TableAppenderConfig(),
    ) as stats_writer:
        for transcript_value in it:
            cdna_seq = transcript_value.transcribe(sequence_func=fv.sequence)
            if len(cdna_seq) == 0:
                continue

            transcript_name = transcript_value.transcript_id
            fa_str = f">{transcript_name}\n{cdna_seq}\n"
            fasta_writer.write(fa_str)
            stats_writer.append(
                (
                    transcript_name,
                    transcript_value.gene_id,
                    transcript_value.seqname,
                    str(transcript_value.start),
                    str(transcript_value.end),
                    transcript_value.strand,
                    str(transcript_value.end - transcript_value.start + 1),
                    str(transcript_value.transcribed_length),
                    str(round(get_gc_percent(cdna_seq) * 100, 2)),
                )
            )
            if write_single_transcript:
                transcript_output_fasta = os.path.join(intermediate_fasta_dir, f"{transcript_name}.fa")
                with get_writer(transcript_output_fasta) as single_transcript_writer:
                    single_transcript_writer.write(fa_str)


def transcribe(
    gt: GeneTreeInterface,
    dst_fasta_path: str,
    fv: FastaViewType,
    show_tqdm: bool = True,
    write_single_transcript: bool = True,
):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    if show_tqdm:
        it = tqdm(iterable=list(gt.transcript_values), desc="Transcribing GTF...")
    else:
        it = gt.transcript_values
    transcribe_transcripts(
        it=it,
        dst_fasta_path=dst_fasta_path,
        fv=fv,
        write_single_transcript=write_single_transcript,
    )


def subset_gtf_by_attribute_value(
    attribute_values: Sequence[str],
    attribute_name: str,
    gtf_filename: str,
    out_filename: str,
    regex: bool = False,
):
    """
    TODO: docs

    .. versionadded:: 1.0.3
    """
    gi = GtfIterator(gtf_filename)
    input_record_num = 0
    intermediate_records = []
    if regex:
        attribute_regex: List[re.Pattern] = list(map(re.compile, attribute_values))
        for gtf_record in gi:
            input_record_num += 1
            this_attribute_value = gtf_record.attribute_get(attribute_name, None)
            if this_attribute_value is None:
                continue
            for possible_regex in attribute_regex:
                if possible_regex.match(this_attribute_value):
                    intermediate_records.append(gtf_record)
                    break
    else:
        attribute_values = list(attribute_values)
        for gtf_record in gi:
            input_record_num += 1
            if gtf_record.attribute_get(attribute_name, None) in attribute_values:
                intermediate_records.append(gtf_record)

    gv = GeneTree.from_feature_iterator(intermediate_records, gene_implementation=DumbGene)
    final_features = list(gv.to_feature_iterator())
    GtfIteratorWriter.write_iterator(final_features, out_filename)
    _lh.info(
        "%d processed with %d (%.2f%%) records output",
        input_record_num,
        len(final_features),
        round(len(final_features) / input_record_num * 100, 2),
    )


def describe(input_filename: str, out_basename: str):
    """
    Describe input GTF.
    """
    gv = GeneTree.from_gtf_file(input_filename, gene_implementation=DumbGene)

    with get_writer(f"{out_basename}.gene.tsv") as gene_writer, get_writer(
        f"{out_basename}.transcripts.tsv"
    ) as transcripts_writer, get_writer(f"{out_basename}.exons.tsv") as exons_writer:
        gene_writer.write(
            "\t".join(
                (
                    "GENE_ID",
                    "TRANSCRIPT_NUMBER",
                    "NAIVE_LENGTH",
                    "TRANSCRIBED_LENGTH",
                    "MAPPABLE_LENGTH",
                    "STRAND",
                )
            )
            + "\n"
        )
        transcripts_writer.write(
            "\t".join(
                (
                    "TRANSCRIPT_ID",
                    "GENE_ID",
                    "NAIVE_LENGTH",
                    "TRANSCRIBED_LENGTH",
                    "EXON_NUMBER",
                    "STRAND",
                )
            )
            + "\n"
        )
        exons_writer.write("\t".join(("TRANSCRIPT_ID", "EXON_NUMBER", "NAIVE_LENGTH", "STRAND")) + "\n")

        for gene in tqdm(
            desc="Iterating over genes...",
            iterable=gv.gene_values,
            total=gv.number_of_genes,
        ):
            gene: Gene
            gene_writer.write(
                "\t".join(
                    (
                        str(gene.gene_id),
                        str(gene.number_of_transcripts),
                        str(gene.naive_length),
                        str(gene.transcribed_length),  # FIXME: Fulfill this requirement
                        str(gene.mappable_length),
                        gene.strand,
                    )
                )
                + "\n"
            )

            transcripts = list(gene.transcript_values)
            for t_i in range(len(transcripts)):
                transcript = transcripts[t_i]
                for exon in list(transcript.exons):
                    exons_writer.write(
                        "\t".join(
                            (
                                exon.transcript_id,
                                str(exon.exon_number),
                                str(exon.naive_length),
                                exon.strand,
                            )
                        )
                        + "\n"
                    )
                transcripts_writer.write(
                    "\t".join(
                        (
                            transcript.transcript_id,
                            transcript.gene_id,
                            str(transcript.naive_length),
                            str(transcript.transcribed_length),
                            str(transcript.number_of_exons),
                            transcript.strand,
                        )
                    )
                    + "\n"
                )
