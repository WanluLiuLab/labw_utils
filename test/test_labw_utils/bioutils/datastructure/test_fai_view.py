import os
import tempfile

import pytest

from labw_utils.bioutils.datastructure.fai_view import FastaIndexView
from labw_utils.bioutils.parser.fai import FastaIndexNotWritableError
from test_labw_utils.bioutils import TEST_DATA_DIR

test_fasta_path = os.path.join(TEST_DATA_DIR, "test.fasta")
test_fai_path = os.path.join(TEST_DATA_DIR, "test.fasta.fai")


def test_fai_io():
    faiv = FastaIndexView.from_fai(test_fai_path)
    assert faiv["chr1"].name == "chr1"
    assert faiv["chr1"].length == 154
    assert faiv["chr1"].offset == 15
    assert faiv["chr1"].line_blen == 27
    assert faiv["chr1"].line_len == 28
    assert len(faiv) == 5

    with tempfile.TemporaryDirectory() as tmpdir:
        out_fq_path = os.path.join(tmpdir, "tmp.fasta.fai")
        faiv.write(out_fq_path)
        assert len(FastaIndexView.from_fai(test_fai_path)) == 5


def test_fasta_to_fai():
    faiv = FastaIndexView.from_fai(test_fai_path)
    faiv_fasta_fh = FastaIndexView.from_fasta(test_fasta_path, full_header=True)
    faiv_fasta_no_fh = FastaIndexView.from_fasta(test_fasta_path, full_header=False)
    assert faiv == faiv_fasta_no_fh
    assert faiv != faiv_fasta_fh
    assert list(faiv.keys()) == ["chr1", "chr2", "chr3", "chr4", "chr6"]
    assert list(faiv_fasta_no_fh.keys()) == ["chr1", "chr2", "chr3", "chr4", "chr6"]
    assert list(faiv_fasta_fh.keys()) == ["chr1 some att", "chr2", "chr3", "chr4", "chr6"]
    with tempfile.TemporaryDirectory() as tmpdir:
        out_fq_path = os.path.join(tmpdir, "tmp.fasta.fai")
        with pytest.raises(FastaIndexNotWritableError):
            faiv_fasta_fh.write(out_fq_path)
        assert not os.path.exists(out_fq_path)
