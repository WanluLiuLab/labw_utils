import os

import numpy as np

from labw_utils.bioutils.record.rmsk_out import parse_record
from test_labw_utils.bioutils import TEST_DATA_DIR

test_rmsk_out_path = os.path.join(TEST_DATA_DIR, "test.rmsk.out")


def test_simple():
    parsed_record = parse_record(
        "  307    0.0  0.0  0.0  chrX             1      262 (17718680) + (CTAAGC)n      Simple_repeat                     1    262     (0)     1  "
    )
    assert parsed_record.sw_score == 307

    parsed_record = parse_record(
        "  267   17.5  1.2 14.1  chrX          8129     8208 (17710734) + PALTTAA2_CE    DNA/PiggyBac                     64    134    (48)    12 *"
    )
    assert parsed_record.sw_score == 267
    assert parsed_record.ssupressed
    parsed_record = parse_record(
        "  410   17.7  0.0 12.8  chrX          5369     5465 (17713477) C PALTTAA2_CE    DNA/PiggyBac                    (4)    178      93     6  "
    )
    assert parsed_record.sw_score == 410
    assert not parsed_record.ssupressed
    assert not parsed_record.strand
    assert parsed_record.sleft == 4
    assert np.all(np.allclose(parsed_record.paln, 0.46703296703296704))
