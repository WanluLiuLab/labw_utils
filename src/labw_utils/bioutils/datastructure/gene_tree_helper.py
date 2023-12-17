"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

import os
import re

from labw_utils.bioutils.algorithm.sequence import get_gc_percent
from labw_utils.bioutils.datastructure.fasta_view import FastaViewType
from labw_utils.bioutils.datastructure.gene_tree import (
    GeneTreeInterface,
    DiploidGeneTree,
)
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.gene import DumbGene, Gene
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.parser.gtf import GtfIterator, GtfIteratorWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.bioutils.record.feature import strand_repr
from labw_utils.commonutils.appender import (
    load_table_appender_class,
    TableAppenderConfig,
    BaseTableAppender,
)
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.file_system import to_safe_filename
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Iterable, Sequence, List, Literal

_lh = get_logger(__name__)


def transcribe_transcripts(
    it: Iterable[Transcript],
    dst_fasta_path: str,
    fv: FastaViewType,
    write_single_transcript: bool = True,
    unsafe_seqname_action: Literal["convert", "error", "skip"] = "convert",
):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """
    intermediate_fasta_dir = ""
    if write_single_transcript:
        intermediate_fasta_dir = dst_fasta_path + ".d"
        os.makedirs(intermediate_fasta_dir, exist_ok=True)
    with FastaWriter(dst_fasta_path) as fasta_writer, load_table_appender_class("TSVTableAppender")(
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
        for transcript in it:
            cdna_seq = transcript.transcribe(sequence_func=fv.sequence)
            if len(cdna_seq) == 0:
                _lh.warning(
                    "seqname '%s' (%s) transcribe failed -- Generated empty cDNA",
                    transcript.transcript_id,
                    repr(transcript),
                )
                continue

            transcript_name = transcript.transcript_id
            safe_transcript_name = to_safe_filename(transcript_name)
            if transcript_name != safe_transcript_name:
                if unsafe_seqname_action == "convert":
                    _lh.warning("seqname '%s' is not safe -- Converted to '%s'", transcript_name, safe_transcript_name)
                elif unsafe_seqname_action == "error":
                    raise ValueError  # TODO
                elif unsafe_seqname_action == "skip":
                    _lh.warning("seqname '%s' is not safe -- skipped", transcript_name)
                    continue

            fasta_writer.write(FastaRecord(seq_id=safe_transcript_name, sequence=cdna_seq))
            stats_writer.append(
                (
                    transcript_name,
                    transcript.gene_id,
                    transcript.seqname,
                    transcript.start,
                    transcript.end,
                    transcript.strand,
                    transcript.naive_length,
                    transcript.transcribed_length,
                    round(get_gc_percent(cdna_seq) * 100, 2),
                )
            )
            if write_single_transcript:
                transcript_output_fasta = os.path.join(intermediate_fasta_dir, f"{transcript_name}.fa")
                with FastaWriter(transcript_output_fasta) as single_transcript_writer:
                    single_transcript_writer.write(FastaRecord(seq_id=safe_transcript_name, sequence=cdna_seq))


def transcribe(
    gt: GeneTreeInterface,
    dst_fasta_path: str,
    fv: FastaViewType,
    show_tqdm: bool = True,
    write_single_transcript: bool = True,
    unsafe_seqname_action: Literal["convert", "error", "skip"] = "convert",
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
        unsafe_seqname_action=unsafe_seqname_action,
    )


def read_partial_gtf_by_attribute_value(
    attribute_values: Sequence[str],
    attribute_name: str,
    gtf_filename: str,
    regex: bool = False,
) -> GeneTreeInterface:
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

    gv = DiploidGeneTree.from_feature_iterator(intermediate_records, gene_implementation=DumbGene)
    final_record_num = len(gv)
    _lh.info(
        "%d processed with %d (%.2f%%) records output",
        input_record_num,
        final_record_num,
        round(final_record_num / input_record_num * 100, 2),
    )
    return gv


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
    gv = read_partial_gtf_by_attribute_value(
        attribute_values=attribute_values,
        attribute_name=attribute_name,
        gtf_filename=gtf_filename,
        regex=regex,
    )
    final_features = list(gv.to_feature_iterator())
    GtfIteratorWriter.write_iterator(final_features, out_filename)


def describe(input_filename: str, out_basename: str):
    """
    Describe input GTF.

    .. versionadd:: 0.1.3
        Migrated from V1API with exon number disabled.
        The file name of the gene description file changed from gene to genes.
    """
    gv = DiploidGeneTree.from_gtf_file(input_filename, gene_implementation=DumbGene)

    with load_table_appender_class("TSVTableAppender")(
        os.path.join(f"{out_basename}.genes"),
        (
            "GENE_ID",
            "TRANSCRIPT_NUMBER",
            "NAIVE_LENGTH",
            "TRANSCRIBED_LENGTH",
            "MAPPABLE_LENGTH",
            "STRAND",
        ),
        TableAppenderConfig(),
    ) as gene_writer, load_table_appender_class("TSVTableAppender")(
        os.path.join(f"{out_basename}.transcripts"),
        (
            "TRANSCRIPT_ID",
            "GENE_ID",
            "NAIVE_LENGTH",
            "TRANSCRIBED_LENGTH",
            "EXON_NUMBER",
            "STRAND",
        ),
        TableAppenderConfig(),
    ) as transcripts_writer, load_table_appender_class(
        "TSVTableAppender"
    )(
        os.path.join(f"{out_basename}.exons"),
        (
            "TRANSCRIPT_ID",
            "EXON_NUMBER",
            "NAIVE_LENGTH",
            "STRAND",
        ),
        TableAppenderConfig(),
    ) as exons_writer:
        gene_writer: BaseTableAppender
        transcripts_writer: BaseTableAppender
        exons_writer: BaseTableAppender
        for gene in tqdm(
            desc="Iterating over genes...",
            iterable=gv.gene_values,
            total=gv.number_of_genes,
        ):
            gene: Gene
            gene_writer.append(
                (
                    gene.gene_id,
                    gene.number_of_transcripts,
                    gene.naive_length,
                    gene.transcribed_length,
                    gene.mappable_length,
                    strand_repr(gene.strand),
                )
            )
            transcripts = list(gene.transcript_values)
            for t_i in range(len(transcripts)):
                transcript = transcripts[t_i]
                for exon in list(transcript.exons):
                    exon: Exon
                    exons_writer.append(
                        (
                            exon.transcript_id,
                            -1,  # This field is disabled
                            exon.naive_length,
                            strand_repr(exon.strand),
                        )
                    )
                transcripts_writer.append(
                    (
                        transcript.transcript_id,
                        transcript.gene_id,
                        transcript.naive_length,
                        transcript.transcribed_length,
                        transcript.number_of_exons,
                        strand_repr(transcript.strand),
                    )
                )
