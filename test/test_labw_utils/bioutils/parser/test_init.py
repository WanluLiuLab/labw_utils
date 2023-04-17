import pytest

from labw_utils.bioutils.parser import FileTypeNotFoundError, BaseFileIterator
from labw_utils.typing_importer import Iterator, Any


class ErrorFileIterator(BaseFileIterator):
    def __iter__(self) -> Iterator[Any]:
        pass


def test_base_file_iterator():
    with pytest.raises(FileTypeNotFoundError):
        _ = ErrorFileIterator("")
