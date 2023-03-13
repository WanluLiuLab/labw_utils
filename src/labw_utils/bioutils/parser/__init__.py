"""
parser -- Basic bioinformatics database parsers

Here contains codes of parsers for basic bioinformatics databases
"""
from abc import abstractmethod, ABC
from typing import Iterable, Optional, IO, Iterator, TypeVar, Generic

from labw_utils.commonutils.stdlib_helper import shutil_helper

_RecordType = TypeVar("_RecordType")


class FileTypeNotFoundError(NameError):
    """
    Developers' Error indicating the subclass was not properly set up.
    """

    def __init__(self, class_name: str):
        super().__init__(f"filetype of class {class_name} was set to None. Please report the bug to developer.")


class _BaseFileIO(ABC):
    _filename: str
    _fd: IO

    filetype: str
    """File type indicator, should be class variable."""

    @property
    def filename(self) -> str:
        """
        Read-only file path.
        """
        return self._filename

    def __new__(cls, *args, **kwargs):
        """
        Method that checks whether a subclass was properly created.
        The subclass should have ``filetype`` attribute.
        """
        _new_instance = super().__new__(cls)
        if not hasattr(_new_instance, "filetype"):
            raise FileTypeNotFoundError(cls.__name__)
        return _new_instance

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


class BaseFileIterator(_BaseFileIO, Generic[_RecordType]):
    """
    Iterate something from a file.
    """

    @abstractmethod
    def __iter__(self) -> Iterator[_RecordType]:
        raise NotImplementedError

    def __init__(self, filename: str, show_tqdm: bool = True, **kwargs):
        _ = kwargs
        self._filename = filename
        self._show_tqdm = show_tqdm


class BaseIteratorWriter(_BaseFileIO, Generic[_RecordType]):

    @staticmethod
    @abstractmethod
    def write_iterator(
            iterable: Iterable[_RecordType],
            filename: str,
            **kwargs
    ):
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
