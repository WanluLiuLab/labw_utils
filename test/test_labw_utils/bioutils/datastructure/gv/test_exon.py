import os

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gv.exon import Exon
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


def test_exon():
    with FastaViewFactory(test_fasta_path, read_into_memory=True) as fasta_view:
        exon = exons[0]
        assert exon.exon_number == 1
        assert exon.transcribe(fasta_view.sequence) == "NNNNNN"
        assert exon.transcript_id == "UN1.1"
        new_exon = exon.update_data(seqname="chr2")
        assert new_exon.seqname == "chr2"
        assert list(new_exon.attribute_keys) == list(exon.attribute_keys)
        assert exons[1] > exons[0]


def test_error_exon():
    with FastaViewFactory(test_fasta_path, read_into_memory=True) as fasta_view:
        exon = Exon(parse_record(
            '1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"'
        ))
        assert exon.exon_number == 0
        assert exon.transcribe(fasta_view.sequence) == ""
        assert exon.transcript_id.startswith("unknown_transcript")
