"""
TODO: docs

.. versionadded:: 1.0.2
"""
from __future__ import annotations

import os

from labw_utils.bioutils.algorithm.sequence import get_gc_percent
from labw_utils.bioutils.datastructure.fasta_view import FastaViewType
from labw_utils.bioutils.datastructure.gene_tree import GeneTreeInterface
from labw_utils.commonutils.appender import load_table_appender_class, TableAppenderConfig
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer


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
        if show_tqdm:
            it = tqdm(iterable=list(gt.transcript_values), desc="Transcribing GTF...")
        else:
            it = gt.transcript_values
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
