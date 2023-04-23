"""
labw_utils.commonutils.lwio -- Enhanced Python IO Functions/Classes

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
TypeError: Type <class '_io.StringIO'> not supported!
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

>>> fio = get_reader(__file__)
>>> fio.__class__.__name__
'TextReaderProxy'
>>> fio.close()
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

import bz2
import enum
import gzip
import io
import lzma
import os
from abc import abstractmethod, ABC

from labw_utils.commonutils.stdlib_helper import shutil_helper
from labw_utils.typing_importer import Iterator, Iterable, List, Union, IO, Optional, AnyStr, Generic, Literal, \
    overload, Dict, Tuple, Type


class SupportsRead(ABC, Generic[AnyStr]):
    @abstractmethod
    def read(self, size: int = 0) -> AnyStr:
        raise NotImplementedError

    def __instancecheck__(self, instance):
        return hasattr(instance, "read")

    def __subclasscheck__(self, subclass):
        return hasattr(subclass, "read")


class SupportsWrite(ABC, Generic[AnyStr]):
    @abstractmethod
    def write(self, size: int = 0) -> AnyStr:
        raise NotImplementedError

    def __instancecheck__(self, instance):
        return hasattr(instance, "write")

    def __subclasscheck__(self, subclass):
        return hasattr(subclass, "write")


PathType = Union[str, bytes, os.PathLike]
FDType = Union[IO, io.IOBase, io.StringIO, io.BytesIO]
PathOrFDType = Union[PathType, FDType]


class ModeEnum(enum.Enum):
    READ = 0
    WRITE = 1
    APPEND = 2


def is_textio(fd: FDType) -> bool:
    """
    Determine whether a given IO is TextIO.

    >>> sio = io.StringIO("")
    >>> is_textio(sio)
    True
    >>> sio.close()

    >>> bio = io.BytesIO(b"")
    >>> is_textio(bio)
    False
    >>> bio.close()

    >>> fio = open(__file__,"r")
    >>> is_textio(fio)
    True
    >>> fio.close()

    >>> fsio = open(__file__,"rt")
    >>> is_textio(fsio)
    True
    >>> fsio.close()

    >>> fbio = open(__file__,"rb")
    >>> is_textio(fbio)
    False
    >>> fbio.close()

    :returns: :py:obj:`True` if is text, :py:obj:`False` (Default) if is binary.
    """
    if isinstance(fd, io.BytesIO):
        return False
    elif isinstance(fd, io.StringIO) or \
            isinstance(fd, io.TextIOBase) or \
            isinstance(fd, io.TextIOWrapper) or \
            hasattr(fd, "encoding"):
        return True
    try:
        modestr = getattr(fd, "mode")
    except AttributeError:
        return False
    if isinstance(modestr, str):
        return "t" in modestr
    return False


def is_io(obj: object) -> bool:
    """
    Determine whether a give object is IO.

    Supports: :py:class:`IO`, :py:class:`io.IOBase` and subclasses.

    >>> fio = open(__file__,"r")
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

    >>> fio = open(__file__,"r")
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


def determine_line_endings(fd: FDType) -> str:
    """
    Determine line endings. If failed, will return OS default.

    This accepts both binary and text IO while return type would always in string.

    :return: One of ``\\r``, ``\\n``, ``\\r\\n``, ``\\n\\r``.
    """
    find_cr = False
    find_lf = False
    while True:
        c = fd.read(1)
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
    _fd: IO[AnyStr]
    """
    The underlying file descriptor.
    """

    def _ensure_open(self):
        """
        Ensure underlying IO is not closed.

        :raises ValueError: If underlying IO is closed.
        """
        if self._fd.closed:
            raise ValueError("I/O operation on closed file.")

    def __instancecheck__(self, instance) -> bool:
        return is_io(instance)

    def __init__(self, fd: FDType):
        """
        Proxy for some file descriptor.

        :raises TypeError: If input cannot be checked with :py:func:`is_io`.
        """
        if not is_io(fd):
            raise TypeError(f"Type {type(fd)} not supported!")
        self._fd = fd  # type: ignore

    @property
    def mode(self) -> str:
        """
        .. warning ::
            This function not always give a string! See following situation:

            >>> import io, gzip
            >>> sio = io.StringIO()
            >>> csio = gzip.open(sio)
            >>> csio.mode
            1
            >>> csio.close()
            >>> sio.close()
        """
        return getattr(self._fd, "mode", "NA")

    @property
    def name(self) -> str:
        return getattr(self._fd, "name", "NA")

    @property
    def closed(self) -> bool:
        return self._fd.closed

    def close(self) -> None:
        self._fd.close()

    def fileno(self) -> int:
        """
        Following are some examples where :py:func:`fileno` are not supported:

        >>> siop = IOProxy(io.StringIO())
        >>> siop.fileno()
        Traceback (most recent call last):
            ...
        TypeError: Type <class '_io.StringIO'> not supported!
        >>> siop.close()

        :return: underlying file descriptor if one exists.
        :raises TypeError: if the IO object does not use a file descriptor.
        """
        try:
            return self._fd.fileno()
        except io.UnsupportedOperation as e:  # In io.StringIO
            raise TypeError(f"Type {type(self._fd)} not supported!") from e

    def flush(self) -> None:
        self._fd.flush()

    def isatty(self) -> bool:
        self._ensure_open()
        return self._fd.isatty()

    def read(self, size: int = -1) -> AnyStr:
        self._ensure_open()
        return self._fd.read(size)

    def readable(self) -> bool:
        self._ensure_open()
        return self._fd.readable()

    def readline(self, limit: int = -1) -> AnyStr:
        self._ensure_open()
        return self._fd.readline(limit)

    def readlines(self, hint: int = -1) -> List[AnyStr]:
        self._ensure_open()
        return self._fd.readlines(hint)

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        self._ensure_open()
        return self._fd.seek(offset, whence)

    def seekable(self) -> bool:
        self._ensure_open()
        return self._fd.seekable()

    def tell(self) -> int:
        self._ensure_open()
        return self._fd.tell()

    def truncate(self, size: Optional[int] = None) -> int:
        self._ensure_open()
        return self._fd.truncate(size)

    def writable(self) -> bool:
        self._ensure_open()
        return self._fd.writable()

    def write(self, s: AnyStr) -> int:
        self._ensure_open()
        return self._fd.write(s)

    def writelines(self, lines: Iterable[AnyStr]) -> None:
        self._ensure_open()
        self._fd.writelines(lines)

    def __next__(self) -> AnyStr:
        self._ensure_open()
        return self._fd.__next__()

    def __iter__(self) -> Iterator[AnyStr]:
        self._ensure_open()
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


class AbstractReader(IOProxy[AnyStr], Generic[AnyStr]):

    def flush(self):
        """Read only"""
        self._ensure_open()

    def readable(self) -> bool:
        self._ensure_open()
        return True

    def truncate(self, *args, **kwargs):
        self._ensure_open()
        raise TypeError("Write operation on Read-Only IOProxy not supported!")

    def writable(self) -> bool:
        self._ensure_open()
        return False

    def write(self, *args, **kwargs):
        self._ensure_open()
        raise TypeError("Write operation on Read-Only IOProxy not supported!")

    def writelines(self, *args, **kwargs):
        self._ensure_open()
        raise TypeError("Write operation on Read-Only IOProxy not supported!")


class TextProxy(IOProxy[str]):

    def __init__(self, fd: FDType):
        super().__init__(fd)
        if not fd.writable():
            raise TypeError(f"Input File Descriptor {fd} not writable!")
        if not fd.readable():
            raise TypeError(f"Input File Descriptor {fd} not readable!")
        if not is_textio(fd):
            raise TypeError(f"Type {type(fd)} is not TextIO!")


class BinaryProxy(IOProxy[bytes]):
    def __init__(self, fd: FDType):
        super().__init__(fd)
        if not fd.writable():
            raise TypeError(f"Input File Descriptor {fd} not writable!")
        if not fd.readable():
            raise TypeError(f"Input File Descriptor {fd} not readable!")
        if is_textio(fd):
            raise TypeError(f"Type {type(fd)} is not BinaryIO!")


class TextReaderProxy(AbstractReader[str]):

    def __init__(self, fd: FDType):
        super().__init__(fd)
        if not fd.readable():
            raise TypeError(f"Input File Descriptor {fd} not readable!")
        if not is_textio(fd):
            raise TypeError(f"Type {type(fd)} is not TextIO!")

    def read(self, size: int = -1) -> str:
        self._ensure_open()
        return self._fd.read(size)

    def readline(self, limit: int = -1) -> str:
        self._ensure_open()
        return self._fd.readline(limit)

    def readlines(self, hint: int = -1) -> List[str]:
        self._ensure_open()
        return self._fd.readlines(hint)

    def __next__(self) -> str:
        self._ensure_open()
        return self._fd.__next__()

    def __iter__(self) -> Iterator[str]:
        self._ensure_open()
        return self._fd.__iter__()


class BinaryReaderProxy(AbstractReader[bytes]):

    def __init__(self, fd: FDType):
        super().__init__(fd)
        if not fd.readable():
            raise TypeError(f"Input File Descriptor {fd} not readable!")
        if is_textio(fd):
            raise TypeError(f"Type {type(fd)} is not BinaryIO!")

    def read(self, size: int = -1) -> bytes:
        self._ensure_open()
        return self._fd.read(size)

    def readline(self, limit: int = -1) -> bytes:
        self._ensure_open()
        return self._fd.readline(limit)

    def readlines(self, hint: int = -1) -> List[bytes]:
        self._ensure_open()
        return self._fd.readlines(hint)

    def __next__(self) -> bytes:
        self._ensure_open()
        return self._fd.__next__()

    def __iter__(self) -> Iterator[bytes]:
        self._ensure_open()
        return self._fd.__iter__()


class AbstractWriter(IOProxy[AnyStr], Generic[AnyStr]):
    def flush(self) -> None:
        self._ensure_open()
        self._fd.flush()

    def read(self, *args, **kwargs):
        self._ensure_open()
        raise TypeError("Read operation on Write-Only IOProxy not supported!")

    def readable(self) -> bool:
        self._ensure_open()
        return False

    def readline(self, *args, **kwargs):
        self._ensure_open()
        raise TypeError("Read operation on Write-Only IOProxy not supported!")

    def readlines(self, *args, **kwargs):
        self._ensure_open()
        raise TypeError("Read operation on Write-Only IOProxy not supported!")

    def writable(self) -> bool:
        self._ensure_open()
        return True

    def write(self, s: AnyStr) -> int:
        self._ensure_open()
        return self._fd.write(s)

    def writelines(self, lines: Iterable[AnyStr]) -> None:
        self._ensure_open()
        self._fd.writelines(lines)

    def __next__(self):
        self._ensure_open()
        raise TypeError("Read operation on Write-Only IOProxy not supported!")

    def __iter__(self):
        self._ensure_open()
        raise TypeError("Read operation on Write-Only IOProxy not supported!")


class BinaryWriterProxy(AbstractWriter[bytes]):

    def __init__(self, fd: FDType):
        super().__init__(fd)
        if not fd.writable():
            raise TypeError(f"Input File Descriptor {fd} not writable!")
        if is_textio(fd):
            raise TypeError(f"Type {type(fd)} is not BinaryIO!")

    def write(self, s: bytes) -> int:
        self._ensure_open()
        return self._fd.write(s)

    def writelines(self, lines: Iterable[bytes]) -> None:
        self._ensure_open()
        self._fd.writelines(lines)


class TextWriterProxy(AbstractWriter[str]):

    def __init__(self, fd: FDType):
        super().__init__(fd)
        if not fd.writable():
            raise TypeError(f"Input File Descriptor {fd} not writable!")
        if not is_textio(fd):
            raise TypeError(f"Type {type(fd)} is not TextIO!")

    def write(self, s: str) -> int:
        self._ensure_open()
        return self._fd.write(s)

    def writelines(self, lines: Iterable[str]) -> None:
        self._ensure_open()
        self._fd.writelines(lines)


class CompressedReaderProxy(AbstractReader):
    @classmethod
    @abstractmethod
    def create(cls, fd: FDType):
        raise NotImplementedError


class CompressedWriterProxy(AbstractWriter):
    @classmethod
    @abstractmethod
    def create(cls, fd: FDType, compression_level: int):
        raise NotImplementedError


class DumbCompressedReaderProxy(CompressedReaderProxy):
    @classmethod
    def create(cls, fd: FDType):
        return cls(fd)


class DumbCompressedWriterProxy(CompressedWriterProxy):
    @classmethod
    def create(cls, fd: FDType, compression_level: int):
        return cls(fd)


class GZipCompressedReaderProxy(CompressedReaderProxy):
    @classmethod
    def create(cls, fd: FDType):
        return cls(gzip.open(fd, "rb"))  # type: ignore


class GZipCompressedWriterProxy(CompressedWriterProxy):
    @classmethod
    def create(cls, fd: FDType, compression_level: int):
        return cls(gzip.open(fd, "wb", compresslevel=compression_level))  # type: ignore


class LZMACompressedReaderProxy(CompressedReaderProxy):
    @classmethod
    def create(cls, fd: FDType):
        return cls(lzma.open(fd, "rb"))  # type: ignore


class LZMACompressedWriterProxy(CompressedWriterProxy):
    @classmethod
    def create(cls, fd: FDType, compression_level: int):
        return cls(lzma.open(fd, "wb", preset=compression_level))  # type: ignore


class BZ2CompressedReaderProxy(CompressedReaderProxy):
    @classmethod
    def create(cls, fd: FDType):
        return cls(bz2.open(fd, "rb"))  # type: ignore


class BZ2CompressedWriterProxy(CompressedWriterProxy):
    @classmethod
    def create(cls, fd: FDType, compression_level: int):
        return cls(bz2.open(fd, "wb", compresslevel=compression_level))  # type: ignore


class CompressionRuleRing:
    _rules: Dict[
        str,
        Tuple[
            Type[CompressedReaderProxy],
            Type[CompressedWriterProxy]
        ]
    ] = {
        "gz": (GZipCompressedReaderProxy, GZipCompressedWriterProxy),
        "gzip": (GZipCompressedReaderProxy, GZipCompressedWriterProxy),
        "GZ": (GZipCompressedReaderProxy, GZipCompressedWriterProxy),
        "lzma": (LZMACompressedReaderProxy, LZMACompressedWriterProxy),
        "xz": (LZMACompressedReaderProxy, LZMACompressedWriterProxy),
        "bz2": (BZ2CompressedReaderProxy, BZ2CompressedWriterProxy),
    }

    @staticmethod
    def get_reader(name: str) -> Type[CompressedReaderProxy]:
        try:
            return CompressionRuleRing._rules[name][0]
        except KeyError:
            return DumbCompressedReaderProxy

    @staticmethod
    def get_writer(name: str) -> Type[CompressedWriterProxy]:
        try:
            return CompressionRuleRing._rules[name][1]
        except KeyError:
            return DumbCompressedWriterProxy


@overload
def file_open(
        file_path: str,
        mode: Literal[ModeEnum.READ],
        is_binary: Literal[False],
        compression: Optional[str] = "Inferred",
        compression_level: int = 0,
        parallel_compression: int = 0
) -> TextReaderProxy:
    ...


@overload
def file_open(
        file_path: str,
        mode: Literal[ModeEnum.WRITE, ModeEnum.APPEND],
        is_binary: Literal[False],
        encoding: str = "UTF-8",
        compression: Optional[str] = "Inferred",
        compression_level: int = 0,
        parallel_compression: int = 0
) -> TextWriterProxy:
    ...


@overload
def file_open(
        file_path: str,
        mode: Literal[ModeEnum.READ],
        is_binary: Literal[True],
        encoding: str = "UTF-8",
        compression: Optional[str] = "Inferred",
        compression_level: int = 0,
        parallel_compression: int = 0
) -> BinaryReaderProxy:
    ...


@overload
def file_open(
        file_path: str,
        mode: Literal[ModeEnum.WRITE, ModeEnum.APPEND],
        is_binary: Literal[True],
        encoding: str = "UTF-8",
        compression: Optional[str] = "Inferred",
        compression_level: int = 0,
        parallel_compression: int = 0
) -> BinaryWriterProxy:
    ...


def file_open(
        file_path: str,
        mode: ModeEnum = ModeEnum.READ,
        is_binary: bool = False,
        encoding: str = "UTF-8",
        compression: Optional[str] = "inferred",
        newline: Optional[str] = None,
        compression_level: int = 0,
        parallel_compression: int = 0
) -> IOProxy:
    if newline is None:
        newline = os.linesep
    if compression == "inferred":
        compression = file_path.split(".")[-1]
    elif compression is None:
        compression = "dumb"
    if mode == ModeEnum.READ:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File {file_path} not found!")
        elif os.path.isdir(file_path):
            raise IsADirectoryError(f"File {file_path} is a directory!")

        rfd = CompressionRuleRing.get_reader(compression).create(
            open(file_path, "rb")
        )
        if is_binary:
            return BinaryReaderProxy(rfd)
        else:
            return TextReaderProxy(io.TextIOWrapper(
                rfd,
                encoding=encoding,
                newline=newline
            ))

    elif mode == ModeEnum.WRITE:
        shutil_helper.touch(file_path)
        wfd = CompressionRuleRing.get_writer(compression).create(
            open(file_path, "wb"),
            compression_level=compression_level
        )
        if is_binary:
            return BinaryWriterProxy(wfd)
        else:
            return TextWriterProxy(io.TextIOWrapper(
                wfd,
                encoding=encoding,
                newline=newline
            ))
    else:
        shutil_helper.touch(file_path)
        afd = CompressionRuleRing.get_writer(compression).create(
            open(file_path, "ab"),
            compression_level=compression_level
        )
        if is_binary:
            return BinaryWriterProxy(afd)
        else:
            return TextWriterProxy(io.TextIOWrapper(
                afd,
                encoding=encoding,
                newline=newline
            ))


def wrap_io(fd: FDType) -> Union[
    TextReaderProxy,
    TextWriterProxy,
    BinaryReaderProxy,
    BinaryWriterProxy,
    TextProxy,
    BinaryProxy
]:
    if is_textio(fd):
        if fd.writable():
            if fd.readable():
                return TextProxy(fd)
            return TextWriterProxy(fd)
        else:
            return TextReaderProxy(fd)
    else:
        if fd.writable():
            if fd.readable():
                return BinaryProxy(fd)
            return BinaryWriterProxy(fd)
        else:
            return BinaryReaderProxy(fd)


@overload
def get_reader(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False],
        **kwargs
) -> IOProxy[bytes]:
    ...


@overload
def get_reader(
        path_or_fd: PathOrFDType,
        is_binary: Literal[True],
        **kwargs
) -> IOProxy[str]:
    ...


@overload
def get_writer(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False],
        **kwargs
) -> IOProxy[bytes]:
    ...


@overload
def get_writer(
        path_or_fd: PathOrFDType,
        is_binary: Literal[True],
        **kwargs
) -> IOProxy[str]:
    ...


@overload
def get_appender(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False],
        **kwargs
) -> IOProxy[bytes]:
    ...


@overload
def get_appender(
        path_or_fd: PathOrFDType,
        is_binary: Literal[True],
        **kwargs
) -> IOProxy[str]:
    ...


def get_reader(
        path_or_fd: PathOrFDType,
        is_binary: bool = False,
        **kwargs
) -> IOProxy:
    if type_check(path_or_fd):
        return wrap_io(path_or_fd)
    else:
        path_or_fd = convert_path_to_str(path_or_fd)
        return file_open(
            path_or_fd,
            mode=ModeEnum.READ,
            is_binary=is_binary,
            **kwargs
        )


def get_writer(
        path_or_fd: PathOrFDType,
        is_binary: bool = False,
        **kwargs
) -> IOProxy:
    if type_check(path_or_fd):
        return wrap_io(path_or_fd)
    else:
        path_or_fd = convert_path_to_str(path_or_fd)
        return file_open(
            path_or_fd,
            mode=ModeEnum.WRITE,
            is_binary=is_binary,
            **kwargs
        )


def get_appender(
        path_or_fd: PathOrFDType,
        is_binary: bool = False,
        **kwargs
) -> IOProxy:
    if type_check(path_or_fd):
        return wrap_io(path_or_fd)
    else:
        path_or_fd = convert_path_to_str(path_or_fd)
        return file_open(
            path_or_fd,
            mode=ModeEnum.APPEND,
            is_binary=is_binary,
            **kwargs
        )
