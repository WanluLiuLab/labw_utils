"""
TODO: docs

.. versionadded:: 1.0.3
"""

__all__ = ("RmskOutIterator",)

from labw_utils.bioutils.parser import BaseFileIterator
from labw_utils.bioutils.record.rmsk_out import RmskOutRecord, parse_record
from labw_utils.commonutils.lwio.safe_io import get_reader
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.typing_importer import Iterator, Final, Iterable


class RmskOutIterator(BaseFileIterator, Iterable[RmskOutRecord]):
    """
    TODO: docs

    .. versionadded:: 1.0.3
    """

    filetype: Final[str] = "RMSK_OUT"
    _skip_first_3_lines: bool

    def __init__(self, filename: str, show_tqdm: bool = True, skip_first_3_lines: bool = True):
        super().__init__(filename, show_tqdm)
        if self._show_tqdm:
            self._fd = get_tqdm_line_reader(self.filename)
        else:
            self._fd = get_reader(self.filename)
        self._skip_first_3_lines = skip_first_3_lines

    def __iter__(self) -> Iterator[RmskOutRecord]:
        if self._skip_first_3_lines:
            for _ in range(3):
                _ = self._fd.readline()
        while True:
            line = self._fd.readline()
            if line == "":
                break

            yield parse_record(line.strip())
        self._fd.close()
