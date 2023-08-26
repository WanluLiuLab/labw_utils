import pytest

from labw_utils.bioutils.record.fastq import (
    FastqRecord,
    SequenceQualityLengthMismatchError,
    MisFormattedFastqRecordError,
)


def test_fastq_errors():
    with pytest.raises(SequenceQualityLengthMismatchError):
        _ = FastqRecord("A", "AA", "AAAA")
    with pytest.raises(MisFormattedFastqRecordError):
        FastqRecord.from_str(["", ""])
    with pytest.raises(MisFormattedFastqRecordError):
        FastqRecord.from_str(["A", "A", "+A", "A"])
    with pytest.raises(MisFormattedFastqRecordError):
        FastqRecord.from_str(["@A", "A", "A", "A"])
