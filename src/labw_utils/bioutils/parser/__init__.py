"""
``labw_utils.bioutils.parser`` -- Basic bioinformatics database parsers

Here contains codes of parsers for basic bioinformatics databases

.. versionadded:: 1.0.2
"""

__all__ = ("BaseFileIterator", "BaseIteratorWriter")

from abc import abstractmethod, ABC

from labw_utils.commonutils.stdlib_helper import shutil_helper
from labw_utils.typing_importer import Iterable, IO, Iterator, TypeVar, Generic

_RecordType = TypeVar("_RecordType")


class _BaseFileIO(ABC):
    _filename: str
    _fd: IO

    @property
    @abstractmethod
    def filetype(self) -> str:
        """File type indicator, should be FINAL class variable."""
        raise NotImplementedError

    @property
    def filename(self) -> str:
        """
        Read-only file path.
        """
        return self._filename

    def __repr__(self) -> str:
        return f"{self.filetype} Iterator for {self._filename} @ {self.tell()}"

    def __str__(self) -> str:
        return repr(self)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.close()

    def tell(self) -> int:
        try:
            return self._fd.tell()
        except (OSError, AttributeError):
            return -1

    def close(self):
        try:
            self._fd.close()
        except AttributeError:
            pass


class BaseFileIterator(_BaseFileIO, Iterable[_RecordType]):
    """
    Iterate something from a file.

    .. versionadded:: 1.0.2
    """

    @abstractmethod
    def __iter__(self) -> Iterator[_RecordType]:
        raise NotImplementedError

    def __init__(self, filename: str, show_tqdm: bool = True, **kwargs):
        _ = kwargs
        self._filename = filename
        self._show_tqdm = show_tqdm


class BaseIteratorWriter(_BaseFileIO, Generic[_RecordType]):
    """
    .. versionadded:: 1.0.2
    """

    @staticmethod
    @abstractmethod
    def write_iterator(iterable: Iterable[_RecordType], filename: str, **kwargs):
        raise NotImplementedError

    def __init__(self, filename: str, **kwargs):
        _ = kwargs
        self._filename = filename

    @abstractmethod
    def write(self, record: _RecordType) -> None:
        raise NotImplementedError

    def destroy_file(self):
        self.close()
        shutil_helper.rm_rf(self._filename)
