import math
import os

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gv.exon import Exon
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.record.gtf import parse_record
from labw_utils.commonutils.stdlib_helper import logger_helper
from test_labw_utils.bioutils import TEST_DATA_DIR
from test_labw_utils.bioutils.datastructure.gv import exons, exon_kwargs

lh = logger_helper.get_logger(__name__)

test_fasta_path = os.path.join(TEST_DATA_DIR, "test.fasta")

standard_transcript = Transcript(
    data=parse_record(
        'chr1\tNA\ttranscript\t5\t35\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.1"'
    ),
    exons=[],
    is_checked=True,
    is_inferred=False,
    keep_sorted=False,
    shortcut=False
)


def test_transcript():
    with FastaViewFactory(test_fasta_path, read_into_memory=True) as fasta_view:
        transcript = Transcript(
            data=standard_transcript.get_data(),
            exons=standard_transcript.exons,
            is_checked=True,
            is_inferred=False,
            keep_sorted=False,
            shortcut=False
        )
        assert transcript.transcript_id == "UN1.1"
        assert transcript.gene_id == "UN1"
        assert transcript.number_of_exons == 0
        assert transcript.transcribed_length == 0
        assert transcript.span_length == -1
        for exon in exons:
            transcript = transcript.add_exon(exon)
        assert transcript.get_exon(0) == exons[0]
        assert transcript.transcribe(fasta_view.sequence) == "NNNNNNNATCGTTACCAT"
        assert transcript.span_length == 26
        assert transcript.naive_length == 31
        assert transcript.exon_boundaries == [(5, 10), (15, 20), (25, 30)]
        assert transcript.splice_sites == [(10, 15), (20, 25)]
        assert transcript.get_intron_length(0) == 5
        assert transcript.get_intron_length(3) == math.inf


def test_infer_names():
    transcript = Transcript(
        data=parse_record(
            'chr1\tNA\ttranscript\t5\t35\t.\t+\t.\tsome_attr "NA";'
        ),
        exons=[],
        is_checked=True,
        is_inferred=False,
        keep_sorted=False,
        shortcut=False
    )
    assert transcript.transcript_id.startswith("unknown_transcript")
    assert transcript.gene_id.startswith("unknown_gene")


def test_transcript_from_exon():
    transcript = Transcript(
        data=parse_record(
            'chr1\tNA\ttranscript\t5\t35\t.\t+\t.\tsome_attr "NA";'
        ),
        exons=[],
        is_checked=True,
        is_inferred=False,
        keep_sorted=False,
        shortcut=False
    )
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
        transcript = transcript.add_exon(exon)
    assert transcript_from_exon.exon_level_equiv(transcript)
    transcript_from_exon = transcript_from_exon.rescale_from_exon_boundaries()
    assert transcript_from_exon.span_length == transcript_from_exon.naive_length


def test_exon_number():
    unnumbered_exons_str = [
        'chr1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.1"; exon_number 0',
        'chr1\tNA\texon\t15\t20\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.2"; exon_number 0',
        'chr1\tNA\texon\t25\t30\t.\t+\t.\tgene_id "UN1"; transcript_id "UN1.3"; exon_number 0'
    ]

    unnumbered_exons = list(map(
        lambda exon_str: Exon(data=parse_record(exon_str), **exon_kwargs),
        unnumbered_exons_str
    ))
    transcript = Transcript(
        data=standard_transcript.get_data(),
        exons=standard_transcript.exons,
        is_checked=True,
        is_inferred=False,
        keep_sorted=False,
        shortcut=False
    )
    for exon in unnumbered_exons:
        transcript = transcript.add_exon(exon)
    updated_transcript = transcript.update_exon_number(exon_number_policy="stranded")
    assert list(exon.attribute_get("exon_number") for exon in transcript.exons) == [0, 0, 0]
    assert list(exon.attribute_get("exon_number") for exon in updated_transcript.exons) == [1, 2, 3]


def test_exon_number_rev():
    unnumbered_exons_str = [
        'chr1\tNA\texon\t5\t10\t.\t-\t.\tgene_id "UN1"; transcript_id "UN1.1"; exon_number 0',
        'chr1\tNA\texon\t15\t20\t.\t-\t.\tgene_id "UN1"; transcript_id "UN1.2"; exon_number 0',
        'chr1\tNA\texon\t25\t30\t.\t-\t.\tgene_id "UN1"; transcript_id "UN1.3"; exon_number 0'
    ]

    unnumbered_exons = list(map(
        lambda exon_str: Exon(data=parse_record(exon_str), **exon_kwargs),
        unnumbered_exons_str
    ))
    transcript = Transcript(
        data=parse_record(
            'chr1\tNA\ttranscript\t5\t35\t.\t-\t.\tgene_id "UN1"; transcript_id "UN1.1"'
        ),
        exons=standard_transcript.exons,
        is_checked=True,
        is_inferred=False,
        keep_sorted=True,
        shortcut=False
    )
    for exon in unnumbered_exons:
        transcript = transcript.add_exon(exon)
    updated_transcript = transcript.update_exon_number(exon_number_policy="stranded")
    assert list(exon.attribute_get("exon_number") for exon in updated_transcript.exons) == [3, 2, 1]
    updated_transcript = transcript.update_exon_number(exon_number_policy="unstranded")
    assert list(exon.attribute_get("exon_number") for exon in updated_transcript.exons) == [1, 2, 3]


def test_del_exon():
    transcript = Transcript(
        data=standard_transcript.get_data(),
        exons=standard_transcript.exons,
        is_checked=True,
        is_inferred=False,
        keep_sorted=True,
        shortcut=False
    )
    for exon in exons:
        transcript = transcript.add_exon(exon)
    transcript_after_del = transcript.del_exon(1)
    assert list(exon.attribute_get("exon_number") for exon in transcript_after_del.exons) == [1, 3]
    transcript_after_del_and_add = transcript_after_del.add_exon(exons[1])
    assert transcript_after_del_and_add.exon_level_equiv(transcript)
    assert not transcript_after_del.exon_level_equiv(transcript)
