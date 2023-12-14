import os
import tempfile

import pandas as pd

from labw_utils.bioutils.datastructure.gene_tree_helper import describe, read_partial_gtf_by_attribute_value
from test_labw_utils.bioutils import TEST_DATA_DIR

test_gtf_path = os.path.join(TEST_DATA_DIR, "test.gtf")


def test_describe():
    with tempfile.TemporaryDirectory() as tmp:
        describe(test_gtf_path, os.path.join(tmp, "tmp"))
        exons = pd.read_csv(os.path.join(tmp, "tmp.exons.tsv"), sep="\t")
        assert len(exons) == 8
        assert list(exons.iloc[0]) == ["UN1.1", -1, 6, "+"]

        transcripts = pd.read_csv(os.path.join(tmp, "tmp.transcripts.tsv"), sep="\t")
        assert len(transcripts) == 4
        assert list(transcripts.iloc[0]) == ["UN1.1", "UN1", 31, 18, 3, "+"]

        genes = pd.read_csv(os.path.join(tmp, "tmp.genes.tsv"), sep="\t")
        assert len(genes) == 3
        assert list(genes.iloc[0]) == ["UN1", 2, 40, 30, 18, "+"]


def test_subset():
    assert (
        read_partial_gtf_by_attribute_value(
            ["UN1.1"],
            "transcript_id",
            test_gtf_path,
            False,
        ).number_of_transcripts
        == 1
    )
    assert (
        read_partial_gtf_by_attribute_value(
            [r"UN1\.*"],
            "transcript_id",
            test_gtf_path,
            True,
        ).number_of_transcripts
        == 2
    )
