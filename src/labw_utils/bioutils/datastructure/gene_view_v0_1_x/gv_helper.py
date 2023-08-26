import os
import re

from labw_utils.bioutils.algorithm.sequence import get_gc_percent
from labw_utils.bioutils.datastructure.fasta_view import FastaViewType
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import (
    GeneViewType,
    GeneViewFactory,
)
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gv_feature_proxy import (
    Transcript,
    Gene,
)
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.old_feature_parser import (
    GtfIterator,
)
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.old_feature_record import (
    GtfRecord,
)
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Tuple, Iterable, Iterator

_lh = get_logger(__name__)


def assert_splice_site_existence(
    this_splice_site: List[Tuple[int, int]],
    all_splice_sites: List[List[Tuple[int, int]]],
) -> bool:
    if this_splice_site in all_splice_sites:
        return True
    else:
        all_splice_sites.append(this_splice_site)
        return False


def get_duplicated_transcript_ids(transcripts: Iterable[Transcript], by_splice_site: bool = True) -> Iterable[str]:
    all_splice_sites: List[List[Tuple[int, int]]] = []

    if by_splice_site:
        for transcript in transcripts:
            this_splice_site = list(transcript.splice_sites)
            if assert_splice_site_existence(this_splice_site, all_splice_sites):
                _lh.warning(f"Will remove {transcript.transcript_id}")
                yield transcript.transcript_id
    else:
        for transcript in transcripts:
            this_splice_site = list(transcript.exon_boundaries)
            if assert_splice_site_existence(this_splice_site, all_splice_sites):
                _lh.warning(f"Will remove {transcript.transcript_id}")
                yield transcript.transcript_id


def gv_dedup(
    gv: GeneViewType,
    by_splice_site: bool = True,
    assume_no_cross_gene_duplication: bool = True,
):
    """
    Remove duplicates in ``transcripts``.

    :param by_splice_site: Detect by splice sites rather than exon boundaries
    :param assume_no_cross_gene_duplication: Whether they may be duplications among genes.
    """
    _lh.info("Finding transcript duplicates in gv...")
    if assume_no_cross_gene_duplication:
        transcript_ids_to_del = []
        for gene in tqdm(iterable=gv.iter_genes()):
            transcript_ids_to_del.extend(
                get_duplicated_transcript_ids(transcripts=gene.iter_transcripts(), by_splice_site=by_splice_site)
            )
    else:
        transcript_ids_to_del = list(
            get_duplicated_transcript_ids(transcripts=gv.iter_transcripts(), by_splice_site=by_splice_site)
        )
    _lh.info(f"Removing {len(transcript_ids_to_del)} transcript duplicate(s) in gv...")
    for transcript_id_to_del in transcript_ids_to_del:
        gv.del_transcript(transcript_id_to_del)
    _lh.info("Removing transcript duplicate(s) in gv FIN")


def transcribe(
    gv: GeneViewType,
    output_fasta: str,
    fv: FastaViewType,
    show_tqdm: bool = True,
    write_single_transcript: bool = True,
):
    intermediate_fasta_dir = output_fasta + ".d"
    os.makedirs(intermediate_fasta_dir, exist_ok=True)
    with get_writer(output_fasta) as fasta_writer, get_writer(output_fasta + ".stats") as stats_writer:
        stats_writer.write(
            "\t".join(
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
                )
            )
            + "\n"
        )
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
            stats_writer.write(
                "\t".join(
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
                + "\n"
            )
            if write_single_transcript:
                transcript_output_fasta = os.path.join(intermediate_fasta_dir, f"{transcript_name}.fa")
                with get_writer(transcript_output_fasta) as single_transcript_writer:
                    single_transcript_writer.write(fa_str)


def subset_gtf_by_attribute_value(
    attribute_values: Iterator[str],
    attribute_name: str,
    gtf_filename: str,
    out_filename: str,
    regex: bool = False,
):
    gi = GtfIterator(gtf_filename)
    input_record_num = 0
    intermediate_records = []
    if regex:
        attribute_regex: List[re.Pattern] = list(map(re.compile, attribute_values))
        for gtf_record in gi:
            input_record_num += 1
            this_attribute_value = gtf_record.attribute.get(attribute_name, None)
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
            if gtf_record.attribute.get(attribute_name, None) in attribute_values:
                intermediate_records.append(gtf_record)

    gv = GeneViewFactory.from_iterable(intermediate_records, record_type=GtfRecord)
    gv.standardize()
    final_record_num = len(gv)
    gv.to_file(out_filename)
    _lh.info(
        "%d processed with %d (%.2f%%) records output",
        input_record_num,
        final_record_num,
        round(final_record_num / input_record_num * 100, 2),
    )


def describe(input_filename: str, out_basename: str):
    """
    Describe input GTF.
    """
    gv = GeneViewFactory.from_file(input_filename, not_save_index=True)

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
            iterable=gv.iter_genes(),
            total=gv.number_of_genes,
        ):
            gene: Gene
            gene_writer.write(
                "\t".join(
                    (
                        str(gene.gene_id),
                        str(gene.number_of_transcripts),
                        str(gene.naive_length),
                        str(gene.transcribed_length),
                        str(gene.mappable_length),
                        gene.strand,
                    )
                )
                + "\n"
            )

            transcripts = list(gene.iter_transcripts())
            for t_i in range(len(transcripts)):
                transcript = transcripts[t_i]
                for exon in list(transcript.iter_exons()):
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
