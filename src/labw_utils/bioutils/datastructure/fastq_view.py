"""
Under construction!
"""
from labw_utils.bioutils.parser.fastq import FastqIterator
from labw_utils.bioutils.record.fastq import FastqRecord
from labw_utils.typing_importer import Dict


class FastqView:
    filename: str
    _dict: Dict[str, FastqRecord]

    def __init__(self, filename: str):
        self.filename = filename
        self._dict = {}
        for fastq_record in FastqIterator(filename):
            self._dict[fastq_record.seq_id] = fastq_record

    def get(self, seq_id: str) -> FastqRecord:
        return self._dict[seq_id]
