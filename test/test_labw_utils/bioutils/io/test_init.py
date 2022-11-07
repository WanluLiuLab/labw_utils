from typing import Iterator, Any

import pytest

from labw_utils.bioutils.io import FileTypeNotFoundError, BaseFileIterator


class ErrorFileIterator(BaseFileIterator):
    def __iter__(self) -> Iterator[Any]:
        pass


def test_base_file_iterator():
    with pytest.raises(FileTypeNotFoundError):
        _ = ErrorFileIterator("")
