import os

from labw_utils.bioutils.record.feature import Feature, FeatureType
from labw_utils.bioutils.record.gtf import parse_record, format_string
from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.typing_importer import List
from test_labw_utils.bioutils import TEST_DATA_DIR

test_gtf_path = os.path.join(TEST_DATA_DIR, "test_various_format.gtf")


def test_gtf_reader():
    with get_reader(test_gtf_path) as reader:
        featl: List[Feature] = list(map(parse_record, reader))
    assert featl[0].seqname == "1"
    assert featl[0].source == "ensembl_havana"
    assert featl[0].feature == "gene"
    assert featl[0].parsed_feature == FeatureType.GENE
    assert featl[0].start == 685679
    assert featl[0].end == 686673
    assert featl[0].score is None
    assert not featl[0].strand
    assert featl[0].frame is None
    assert list(featl[0].attribute_keys) == ["gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype"]
    assert featl[0].attribute_get("gene_version") == 1
    assert featl[0].attribute_get("not_exist") is None

    assert featl[1].score == 1024.3
    assert featl[1].strand is None
    assert featl[1].parsed_feature is FeatureType.NOT_PRESENT
    assert list(featl[1].attribute_values) == ["DDX11L1", "NR_046018", "DDX11L1"]

    with get_reader(test_gtf_path) as reader:
        assert list(map(format_string, featl)) == list(map(str.strip, reader.readlines()))


def test_partial_gtf_reader():
    with get_reader(test_gtf_path) as reader:
        featl: List[Feature] = list(map(
            lambda gtf_str: parse_record(
                gtf_str,
                skip_fields=["score", "not_exist"],
                included_attributes=["gene_name", "not_exist"]
            ),
            reader
        ))
    assert featl[0].score is None
    assert list(featl[1].attribute_keys) == ["gene_name"]


def test_format_str():
    gtf_str = '1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"; exon_number 1; description "Sample exon"'
    assert format_string(parse_record(gtf_str), "all") == \
           '1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"; exon_number "1"; description "Sample exon";'
    assert format_string(parse_record(gtf_str), "none") == \
           '1\tNA\texon\t5\t10\t.\t+\t.\tgene_id UN1; exon_number 1; description Sample exon;'
    assert format_string(parse_record(gtf_str), "blank") == \
           '1\tNA\texon\t5\t10\t.\t+\t.\tgene_id UN1; exon_number 1; description "Sample exon";'
    assert format_string(parse_record(gtf_str), "string") == \
           '1\tNA\texon\t5\t10\t.\t+\t.\tgene_id "UN1"; exon_number 1; description "Sample exon";'
