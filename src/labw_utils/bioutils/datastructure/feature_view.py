"""
Under construction!
"""

from typing import Optional, List

from labw_utils.bioutils.parser.feature import GtfIterator
from labw_utils.bioutils.record.feature import GtfRecord


class GtfView(List[GtfRecord]):
    filename: str

    def __init__(self, filename: Optional[str] = None):
        super().__init__()
        if filename is None:
            self.filename = "in_memory"
        else:
            self.filename = filename
            self.load()

    def load(self):
        self.clear()
        for gtf_record in GtfIterator(self.filename):
            self.append(gtf_record)
