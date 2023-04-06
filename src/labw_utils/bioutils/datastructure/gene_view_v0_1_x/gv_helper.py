import os
from typing import List, Tuple, Iterable

from labw_utils.bioutils.algorithm.sequence import get_gc_percent
from labw_utils.bioutils.datastructure.fasta_view import FastaViewType
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewType
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_feature_proxy import Transcript
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def assert_splice_site_existence(
        this_splice_site: List[Tuple[int, int]],
        all_splice_sites: List[List[Tuple[int, int]]]
) -> bool:
    if this_splice_site in all_splice_sites:
        return True
    else:
        all_splice_sites.append(this_splice_site)
        return False


def get_duplicated_transcript_ids(
        transcripts: Iterable[Transcript],
        by_splice_site: bool = True
) -> Iterable[str]:
    all_splice_sites: List[List[Tuple[int, int]]] = []

    if by_splice_site:
        for transcript in transcripts:
            this_splice_site = list(transcript.splice_sites)
            if assert_splice_site_existence(this_splice_site, all_splice_sites):
                lh.warning(f"Will remove {transcript.transcript_id}")
                yield transcript.transcript_id
    else:
        for transcript in transcripts:
            this_splice_site = list(transcript.exon_boundaries)
            if assert_splice_site_existence(this_splice_site, all_splice_sites):
                lh.warning(f"Will remove {transcript.transcript_id}")
                yield transcript.transcript_id


def gv_dedup(
        gv: GeneViewType,
        by_splice_site: bool = True,
        assume_no_cross_gene_duplication: bool = True
):
    """
    Remove duplicates in ``transcripts``.

    :param by_splice_site: Detect by splice sites rather than exon boundaries
    :param assume_no_cross_gene_duplication: Whether they may be duplications among genes.
    """
    lh.info("Finding transcript duplicates in gv...")
    if assume_no_cross_gene_duplication:
        transcript_ids_to_del = []
        for gene in tqdm(iterable=gv.iter_genes()):
            transcript_ids_to_del.extend(get_duplicated_transcript_ids(
                transcripts=gene.iter_transcripts(),
                by_splice_site=by_splice_site
            ))
    else:
        transcript_ids_to_del = list(get_duplicated_transcript_ids(
            transcripts=gv.iter_transcripts(),
            by_splice_site=by_splice_site
        ))
    lh.info(f"Removing {len(transcript_ids_to_del)} transcript duplicate(s) in gv...")
    for transcript_id_to_del in transcript_ids_to_del:
        gv.del_transcript(transcript_id_to_del)
    lh.info("Removing transcript duplicate(s) in gv FIN")


def transcribe(
        gv: GeneViewType,
        output_fasta: str,
        fv: FastaViewType,
        show_tqdm: bool = True
):
    intermediate_fasta_dir = output_fasta + ".d"
    os.makedirs(intermediate_fasta_dir, exist_ok=True)
    with get_writer(output_fasta) as fasta_writer, \
            get_writer(output_fasta + ".stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "GENE_ID",
            "SEQNAME",
            "START",
            "END",
            "STRAND",
            "ABSOLUTE_LENGTH",
            "TRANSCRIBED_LENGTH",
            "GC"
        )) + "\n")
        if show_tqdm:
            it = tqdm(iterable=list(gv.iter_transcripts()), desc="Transcribing GTF...")
        else:
            it = gv.iter_transcripts()
        for transcript_value in it:
            transcript_value: Transcript
            cdna_seq = transcript_value.cdna_sequence(sequence_func=fv.sequence)
            if len(cdna_seq) == 0:
                continue

            transcript_name = transcript_value.transcript_id
            fa_str = f">{transcript_name}\n{cdna_seq}\n"
            fasta_writer.write(fa_str)
            stats_writer.write("\t".join((
                transcript_name,
                transcript_value.gene_id,
                transcript_value.seqname,
                str(transcript_value.start),
                str(transcript_value.end),
                transcript_value.strand,
                str(transcript_value.end - transcript_value.start + 1),
                str(transcript_value.transcribed_length),
                str(round(get_gc_percent(cdna_seq) * 100, 2))
            )) + "\n")
            transcript_output_fasta = os.path.join(intermediate_fasta_dir, f"{transcript_name}.fa")
            with get_writer(transcript_output_fasta) as single_transcript_writer:
                single_transcript_writer.write(fa_str)
