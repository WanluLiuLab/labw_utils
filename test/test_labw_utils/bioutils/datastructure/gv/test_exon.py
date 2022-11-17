import os

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.record.gtf import parse_record
from labw_utils.commonutils.stdlib_helper import logger_helper
from test_labw_utils.bioutils import TEST_DATA_DIR
from test_labw_utils.bioutils.datastructure.gv import exons

lh = logger_helper.get_logger(__name__)

test_fasta_path = os.path.join(TEST_DATA_DIR, "test.fasta")


def test_exon():
    with FastaViewFactory(test_fasta_path, read_into_memory=True) as fasta_view:
        exon = exons[0]
        assert exon.exon_number == 1
        assert exon.transcribe(fasta_view.sequence) == "NNNNNN"
        assert len(exon.transcribe(fasta_view.sequence)) == exon.transcribed_length
        assert exon.transcribed_length == exon.naive_length
        assert exon.transcript_id == "UN1.1"
        new_exon = Exon(data=exon.get_data().update(seqname="chr2"), is_checked=False, shortcut=True)
        assert new_exon.seqname == "chr2"
        assert list(new_exon.attribute_keys) == list(exon.attribute_keys)
        assert exons[1] > exons[0]


def test_error_exon():
    with FastaViewFactory(test_fasta_path, read_into_memory=True) as fasta_view:
        exon = Exon(
            data=parse_record(
                '1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"'
            ),
            is_checked=True,
            shortcut=False
        )
        assert exon.exon_number == 0
        assert exon.transcribe(fasta_view.sequence) == ""
        assert exon.transcript_id.startswith("unknown_transcript")
        assert exon.gene_id == "UN1"
