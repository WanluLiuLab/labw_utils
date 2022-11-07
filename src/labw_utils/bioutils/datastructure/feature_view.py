"""
Under construction!
"""

from typing import Optional, List

from labw_utils.bioutils.io.feature import GtfIterator
from labw_utils.bioutils.typing.feature import GtfRecord


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
