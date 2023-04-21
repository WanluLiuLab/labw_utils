"""
labw_utils.commonutils.io -- Enhanced Python IO Functions/Classes

This module includes some enhanced IO functions,
like IO that automatically creates missing intermediate directories or IO with tqdm progress bar.

It also supports Python standard archive types like ``bz2``, ``gzip`` or ``lzma``
like is implemented in :py:mod:`fileinput` module.

Following is an example using :py:class:`IOProxy` over :py:class:`io.StringIO`.

>>> sio = io.StringIO("AAAA")
>>> siop = IOProxy(sio)
>>> siop.mode
'NA'
>>> siop.name
'NA'
>>> siop.closed
False
>>> siop.fileno()
Traceback (most recent call last):
    ...
OSError
>>> siop.isatty()
False
>>> siop.seek(0)
0
>>> siop.tell()
0
>>> siop.readlines()
['AAAA']
>>> siop.writable()
True
>>> siop.seek(0)
0
>>> list(siop)
['AAAA']
>>> siop.close()
>>> siop.closed
True
>>> sio.closed
True
"""

from __future__ import annotations

__all__ = (
    "determine_line_endings",
    "PathType",
    "PathOrFDType",
    "FDType",
    "IOProxy",
    "is_io",
    "is_path",
    "type_check"
)

import enum
import io
import os

from labw_utils.typing_importer import Iterator, Iterable, List, Union, IO, Optional, AnyStr, Generic, Literal, overload

PathType = Union[str, bytes, os.PathLike]
FDType = Union[IO, io.IOBase, io.StringIO, io.BytesIO]
PathOrFDType = Union[PathType, FDType]
class ModeEnum(enum.Enum):
    READ = 0
    WRITE = 1
    APPEND = 2


def is_io(obj: object) -> bool:
    """
    Determine whether a give object is IO.

    Supports: :py:class:`IO`, :py:class:`io.IOBase` and subclasses.

    >>> fio = open(__file__, "r")
    >>> is_io(fio)
    True
    >>> fio.close()

    >>> sio = io.StringIO("AAAA")
    >>> is_io(sio)
    True
    >>> sio.close()

    >>> bio = io.BytesIO(b"AAAA")
    >>> is_io(bio)
    True
    >>> bio.close()

    >>> is_io("AAA")
    False
    """
    if (
            isinstance(obj, IO) or
            isinstance(obj, io.IOBase)
    ):
        return True
    return False


def is_path(obj: object) -> bool:
    """
    Check whether a given object is path.

    Supports: :py:class:`str`, :py:class:`bytes`, :py:class:`os.PathLike` and subclasses.

    >>> import pathlib
    >>> is_path("/")
    True
    >>> is_path(b"/")
    True
    >>> is_path(pathlib.Path())
    True
    >>> is_path(True)
    False
    """
    if (
            isinstance(obj, os.PathLike) or
            isinstance(obj, bytes) or
            isinstance(obj, str)
    ):
        return True
    return False


def convert_path_to_str(path: PathType) -> str:
    """
    Convert path to string.

    >>> import pathlib
    >>> convert_path_to_str("/")
    '/'
    >>> convert_path_to_str(b"/")
    '/'
    >>> convert_path_to_str(pathlib.Path())
    '.'

    :raises TypeError: On unexpected types.
    """
    if isinstance(path, os.PathLike):
        return path.__fspath__()
    elif isinstance(path, bytes):
        return str(path, encoding="UTF-8")
    elif isinstance(path, str):
        return path
    else:
        raise TypeError(f"Type {type(path)} not supported!")


def type_check(obj: object) -> bool:
    """
    Check whether input object is path or file descriptor.

    >>> import pathlib

    >>> fio = open(__file__, "r")
    >>> type_check(fio)
    True
    >>> fio.close()

    >>> sio = io.StringIO("AAAA")
    >>> type_check(sio)
    True
    >>> sio.close()

    >>> bio = io.BytesIO(b"AAAA")
    >>> type_check(bio)
    True
    >>> bio.close()

    >>> type_check("/")
    False
    >>> type_check(b"/")
    False
    >>> type_check(pathlib.Path())
    False
    >>> type_check(True)
    Traceback (most recent call last):
        ...
    TypeError: Type <class 'bool'> not supported!

    :return: :py:obj:`True` if is IO, :py:obj:`False` if is string.
    :raises TypeError: On unexpected types.
    """
    if is_io(obj):
        return True
    elif is_path(obj):
        return False
    else:
        raise TypeError(f"Type {type(obj)} not supported!")


def determine_line_endings(file_path_or_fd: PathOrFDType) -> str:
    """
    Determine line endings. If failed, will return OS default.

    This accepts both binary and text IO while return type would always in string.

    :param file_path_or_fd: File path or file descriptor.
    :return: One of ``\\r``, ``\\n``, ``\\r\\n``, ``\\n\\r``.
    """
    if not type_check(file_path_or_fd):
        file_path_or_fd = convert_path_to_str(file_path_or_fd)
        return determine_line_endings(open(file_path_or_fd, "r"))

    find_cr = False
    find_lf = False
    while True:
        c = file_path_or_fd.read(1)
        if c is None:
            break
        elif c == '\r' or c == b'\r':
            find_cr = True
            if find_lf:
                return '\n\r'
        elif c == '\n' or c == b'\n':
            find_lf = True
            if find_cr:
                return '\r\n'
        else:
            if find_cr:
                return '\r'
            if find_lf:
                return '\n'
            find_cr = False
            find_lf = False
    return os.linesep


class IOProxy(IO[AnyStr]):
    """
    IO Proxy for IO objects.
    """
    _fd: FDType
    """
    The underlying file descriptor.
    """

    def __instancecheck__(self, instance) -> bool:
        return is_io(instance)

    def __init__(self, fd: FDType, *args, **kwargs):
        """
        Proxy for some file descriptor.

        :raises TypeError: If input cannot be checked with :py:func:`is_io`.
        """
        if not is_io(fd):
            raise TypeError(f"Type {type(fd)} not supported!")
        self._fd = fd

    @property
    def mode(self) -> str:
        if hasattr(self._fd, 'mode'):
            return self._fd.mode
        else:
            return "NA"

    @property
    def name(self) -> str:
        if hasattr(self._fd, 'name'):
            return self._fd.name
        else:
            return "NA"

    @property
    def closed(self) -> bool:
        return self._fd.closed

    def close(self) -> None:
        self._fd.close()

    def fileno(self) -> int:
        """
        :return: underlying file descriptor if one exists.
        :raises OSError: if the IO object does not use a file descriptor.
        """
        try:
            return self._fd.fileno()
        except io.UnsupportedOperation as e:  # In io.StringIO
            raise OSError from e

    def flush(self) -> None:
        self._fd.flush()

    def isatty(self) -> bool:
        return self._fd.isatty()

    def read(self, size: int = -1) -> AnyStr:
        return self._fd.read(size)

    def readable(self) -> bool:
        return self._fd.readable()

    def readline(self, limit: int = -1) -> AnyStr:
        return self._fd.readline(limit)

    def readlines(self, hint: int = -1) -> List[AnyStr]:
        return self._fd.readlines(hint)

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        return self._fd.seek(offset, whence)

    def seekable(self) -> bool:
        return self._fd.seekable()

    def tell(self) -> int:
        return self._fd.tell()

    def truncate(self, size: Optional[int] = None) -> int:
        return self._fd.truncate(size)

    def writable(self) -> bool:
        return self._fd.writable()

    def write(self, s: AnyStr) -> int:
        return self._fd.write(s)

    def writelines(self, lines: Iterable[AnyStr]) -> None:
        self._fd.writelines(lines)

    def __next__(self) -> AnyStr:
        return self._fd.__next__()

    def __iter__(self) -> Iterator[AnyStr]:
        return self._fd.__iter__()

    def __enter__(self):
        try:
            self._fd.__enter__()
        except AttributeError:
            pass
        return self

    def __exit__(self, *args, **kwargs):
        try:
            self._fd.__exit__(*args, **kwargs)
        except AttributeError:
            return
