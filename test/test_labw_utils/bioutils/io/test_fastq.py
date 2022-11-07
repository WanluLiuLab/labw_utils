import os
import random
import tempfile
from typing import List

from labw_utils.bioutils.io.fastq import FastqIterator, FastqWriter
from labw_utils.bioutils.typing.fastq import FastqRecord
from test_labw_utils.bioutils import TEST_DATA_DIR

test_fastq_path = os.path.join(TEST_DATA_DIR, "test.fastq")


def test():
    fql: List[FastqRecord] = list(iter(FastqIterator(test_fastq_path)))
    assert fql[0].seq_id == "NR_138525.1_318_626_0_1_0_0_1:0:0_3:0:0_0/2"
    assert fql[0].sequence == "TTCAAATAAAGCACTTTATTGACCAGTTAAAATATCCAGAAGAGTTGACA"
    assert fql[0].quality == "5321421222012227124162/2223/2125302220312213224522"
    assert len(fql[0]) == 50
    assert len(fql) == 5

    sampled_fql = random.sample(fql, 3)
    with tempfile.TemporaryDirectory() as tmpdir:
        out_fq_path = os.path.join(tmpdir, "tmp.fq")
        FastqWriter.write_iterator(sampled_fql, out_fq_path)
        assert len(list(iter(FastqIterator(out_fq_path)))) == 3
