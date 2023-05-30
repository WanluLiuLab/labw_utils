import os
import random
import tempfile

import pytest

from labw_utils.bioutils.parser.fai import FastaBasedFastaIndexIterator, FAIBasedFastaIndexIterator, \
    FastaIndexNotWritableError, FastaIndexWriter, DuplicatedFastaNameError
from labw_utils.bioutils.record.fai import FastaIndexRecord
from labw_utils.typing_importer import List
from test_labw_utils import NULL_PATH
from test_labw_utils.bioutils import TEST_DATA_DIR

test_fasta_path = os.path.join(TEST_DATA_DIR, "test.fasta")
test_dup_fasta_path = os.path.join(TEST_DATA_DIR, "test_duplicate.fasta")
test_fai_path = os.path.join(TEST_DATA_DIR, "test.fasta.fai")


def test_fai_io():
    fairl: List[FastaIndexRecord] = list(iter(FAIBasedFastaIndexIterator(test_fai_path)))
    assert fairl[0].name == 'chr1'
    assert fairl[0].length == 154
    assert fairl[0].offset == 15
    assert fairl[0].line_blen == 27
    assert fairl[0].line_len == 28
    assert len(fairl) == 5

    sampled_fql = random.sample(fairl, 3)
    with tempfile.TemporaryDirectory() as tmpdir:
        out_fq_path = os.path.join(tmpdir, "tmp.fasta.fai")
        FastaIndexWriter.write_iterator(sampled_fql, out_fq_path)
        assert len(list(iter(FAIBasedFastaIndexIterator(out_fq_path)))) == 3


def test_fasta_to_fai():
    fairl: List[FastaIndexRecord] = list(iter(FAIBasedFastaIndexIterator(test_fai_path)))
    fairl_fasta_fh: List[FastaIndexRecord] = list(
        iter(FastaBasedFastaIndexIterator(test_fasta_path, full_header=True))
    )
    fairl_fasta_no_fh: List[FastaIndexRecord] = list(
        iter(FastaBasedFastaIndexIterator(test_fasta_path, full_header=False))
    )
    assert fairl == fairl_fasta_no_fh
    assert fairl != fairl_fasta_fh
    with tempfile.TemporaryDirectory() as tmpdir:
        out_fq_path = os.path.join(tmpdir, "tmp.fasta.fai")
        with pytest.raises(FastaIndexNotWritableError):
            FastaIndexWriter.write_iterator(fairl_fasta_fh, out_fq_path)
        assert not os.path.exists(out_fq_path)


def test_duplicated_fasta():
    with pytest.raises(DuplicatedFastaNameError):
        _ = list(FastaBasedFastaIndexIterator(test_dup_fasta_path, full_header=False))


def test_empty_file():
    assert len(list(FastaBasedFastaIndexIterator(NULL_PATH))) == 0
    assert len(list(FAIBasedFastaIndexIterator(NULL_PATH))) == 0
