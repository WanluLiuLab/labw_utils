import os

from labw_utils.bioutils.datastructure.gene_tree import GeneTree
from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.commonutils.stdlib_helper import logger_helper
from test_labw_utils.bioutils import TEST_DATA_DIR

lh = logger_helper.get_logger(__name__)

test_gtf_path = os.path.join(TEST_DATA_DIR, "test.gtf")


# FIXME
def test():
    gt = GeneTree.from_feature_iterator(GtfIterator(test_gtf_path))
    assert list(gt.transcript_ids) == ["UN1.1", "UN1.2", "UN2.1", "UN3.1"]
    assert list(gt.gene_ids) == ["UN1", "UN2", "UN3"]
    assert list(gt.get_gene("UN1")[0].transcript_ids) == ["UN1.1", "UN1.2"]
    assert list(gt.get_gene("UN2")[0].transcript_ids) == ["UN2.1"]
    assert list(gt.get_gene("UN3")[0].transcript_ids) == ["UN3.1"]
    gt = gt.replace_transcript(gt.get_transcript("UN2.1").update_exon_number(exon_number_policy="stranded"))
    assert [exon.attribute_get("exon_number") for exon in gt.get_transcript("UN2.1").exons] == [2, 1]
    gt = gt.replace_transcript(gt.get_transcript("UN2.1").update_exon_number(exon_number_policy="unstranded"))
    assert [exon.attribute_get("exon_number") for exon in gt.get_transcript("UN2.1").exons] == [1, 2]


if __name__ == "__main__":
    test()
