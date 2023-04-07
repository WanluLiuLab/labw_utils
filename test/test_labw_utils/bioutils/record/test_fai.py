import pytest

from labw_utils.bioutils.record.fai import FastaIndexRecord, FastaIndexRecordParserError


def test():
    fai_str = "chr1\t154\t15\t27\t28"
    fair = FastaIndexRecord.from_fai_str(fai_str)
    assert fair.name == 'chr1'
    assert fair.length == 154
    assert fair.offset == 15
    assert fair.line_blen == 27
    assert fair.line_len == 28
    assert str(fair) == fai_str

    with pytest.raises(FastaIndexRecordParserError):
        _ = FastaIndexRecord.from_fai_str("AAAAA")
