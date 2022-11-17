import os

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.record.gtf import parse_record
from labw_utils.commonutils.stdlib_helper import logger_helper
from test_labw_utils.bioutils import TEST_DATA_DIR
from test_labw_utils.bioutils.datastructure.gv import exons

lh = logger_helper.get_logger(__name__)

test_fasta_path = os.path.join(TEST_DATA_DIR, "test.fasta")


def test_transcript():
    with FastaViewFactory(test_fasta_path, read_into_memory=True) as fasta_view:
        transcript = Transcript(
            data=parse_record(
                'chr1\tNA\ttranscript\t5\t35\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.1"'
            ),
            exons=[],
            is_checked=True,
            is_inferred=False,
            keep_sorted=False,
            shortcut=False
        )
        for exon in exons:
            transcript = transcript.add_exon(exon)
        assert transcript.transcribe(fasta_view.sequence) == "NNNNNNNATCGTTACCAT"
        assert transcript.span_length == 26
        assert transcript.naive_length == 31

        transcript_from_exon = Transcript(
            data=exons[0].get_data(),
            exons=[],
            is_checked=True,
            is_inferred=True,
            keep_sorted=False,
            shortcut=False
        )
        for exon in exons:
            transcript_from_exon = transcript_from_exon.add_exon(exon)
        assert transcript_from_exon.exon_level_equiv(transcript)
        transcript_from_exon = transcript_from_exon.rescale_from_exon_boundaries()
        assert transcript_from_exon.span_length == transcript_from_exon.naive_length
