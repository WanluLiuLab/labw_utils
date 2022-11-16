import os

import pytest

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.transcript import ExonInATranscriptOnDifferentChromosomeError, \
    DuplicatedExonError, ExonInATranscriptOnDifferentStrandError, Transcript
from labw_utils.bioutils.record.gtf import parse_record
from labw_utils.commonutils.stdlib_helper import logger_helper
from test_labw_utils.bioutils import TEST_DATA_DIR

lh = logger_helper.get_logger(__name__)

test_fasta_path = os.path.join(TEST_DATA_DIR, "test.fasta")

exons = [
    Exon(parse_record(
        'chr1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.1"; exon_number 1'
    )),
    Exon(parse_record(
        'chr1\tNA\texon\t15\t20\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.2"; exon_number 2'
    )),
    Exon(parse_record(
        'chr1\tNA\texon\t25\t30\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.3"; exon_number3'
    ))
]


def test_transcript():
    with FastaViewFactory(test_fasta_path, read_into_memory=True) as fasta_view:
        transcript = Transcript(parse_record(
            'chr1\tNA\ttranscript\t5\t35\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.1"'
        ))
        for exon in exons:
            transcript = transcript.add_exon(exon)
        assert transcript.transcribe(fasta_view.sequence) == "NNNNNNNATCGTTACCAT"
        with pytest.raises(DuplicatedExonError):
            transcript.add_exon(exons[0])
        with pytest.raises(DuplicatedExonError):
            transcript.add_exon(exons[1])
        with pytest.raises(DuplicatedExonError):
            transcript.add_exon(exons[2])
        with pytest.raises(ExonInATranscriptOnDifferentChromosomeError):
            transcript.add_exon(exons[2].update_data(seqname="chr2"))
        with pytest.raises(ExonInATranscriptOnDifferentStrandError):
            transcript.add_exon(exons[2].update_data(strand="."))
        assert transcript.span_length == 26
        assert transcript.naive_length == 31

        transcript_from_exon = Transcript.infer_from_exon(exons[0])
        for exon in exons:
            transcript_from_exon = transcript_from_exon.add_exon(exon)
        assert transcript_from_exon.exon_level_equiv(transcript)
        transcript_from_exon = transcript_from_exon.rescale_from_exon_boundaries()
        assert transcript_from_exon.span_length == transcript_from_exon.naive_length
