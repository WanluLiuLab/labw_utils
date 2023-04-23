import os

from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.feature import FeatureType
from labw_utils.bioutils.record.gtf import format_string
from labw_utils.commonutils.lwio.safe_io import get_reader
from test_labw_utils.bioutils import TEST_DATA_DIR

test_gtf_path = os.path.join(TEST_DATA_DIR, "test_various_format.gtf")


def test_gtf_reader():
    with GtfIterator(test_gtf_path) as gtfi:
        featl = list(gtfi)
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

# TODO: Writer test
