import pytest

from labw_utils.bioutils.parser import BaseFileIterator
from labw_utils.typing_importer import Iterator, Any


class ErrorFileIterator(BaseFileIterator):

    @property
    def filetype(self) -> str:
        return ""

    def __iter__(self) -> Iterator[Any]:
        pass
