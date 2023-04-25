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

import bz2
import enum
import gzip
import io
import lzma
import os
from abc import abstractmethod

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.typing_importer import Iterator, Iterable, List, Union, IO, Optional, AnyStr, Generic, Literal, \
    overload, Dict, Tuple, Type


def wc_l_io(fd: FDType, block_size: int = 4096) -> int:
    """
    Count lines in a file.

    :param fd: A finite seekable block IO object.
    :param block_size: Number of bytes to read at once.
    :return: Line number. -1 if not seekable.
    """
    if fd.seekable():
        curr_pos = fd.tell()
    else:
        return -1
    reti = 0
    fd.seek(0)
    block = fd.read(block_size)
    """Assume 4k aligned filesystem"""
    if len(block) == 0:
        return 0
    else:
        if isinstance(block, str):
            reti += block.count("\n")
            while True:
                block = fd.read(block_size)
                if len(block) == 0:
                    break
                reti += block.count("\n")  # type: ignore
        else:
            reti += block.count(b"\n")
            while True:
                block = fd.read(block_size)
                if len(block) == 0:
                    break
                reti += block.count(b"\n")  # type: ignore
    if fd.tell() != 0 and reti == 0:
        reti = 1  # To keep similar behaviour to GNU WC
    fd.seek(curr_pos)
    return reti


def wc_c_io(fd: FDType) -> int:
    """
    Count the number of chars inside a file, i.e. File length.

    :param fd: A finite seekable IO object.
    :return: File length. -1 if not seekable.
    """
    if fd.seekable():
        curr_pos = fd.tell()
        fd.seek(0, 2)
        reti = fd.tell()
        fd.seek(curr_pos)
        return reti
    else:
        return -1


# class SupportsRead(ABC, Generic[AnyStr]):
#     @abstractmethod
#     def read(self, size: int = 0) -> AnyStr:
#         raise NotImplementedError
#
#     def __instancecheck__(self, instance):
#         return hasattr(instance, "read")
#
#     def __subclasscheck__(self, subclass):
#         return hasattr(subclass, "read")
#
#
# class SupportsWrite(ABC, Generic[AnyStr]):
#     @abstractmethod
#     def write(self, size: int = 0) -> AnyStr:
#         raise NotImplementedError
#
#     def __instancecheck__(self, instance):
#         return hasattr(instance, "write")
#
#     def __subclasscheck__(self, subclass):
#         return hasattr(subclass, "write")


PathType = Union[str, bytes, os.PathLike]
FDType = Union[IO, io.IOBase, io.StringIO, io.BytesIO]
PathOrFDType = Union[PathType, FDType]


class ModeEnum(enum.Enum):
    """Simple read- or write-only mode enum."""
    READ = 0
    """To open a file read-only"""
    WRITE = 1
    """To open a file write-only. Will create a new file if not exist."""
    APPEND = 2
    """To open a file write-only without truncating contents inside. Will create a new file if not exist."""


def io_repr(fd: FDType):
    """
    >>> sio = io.StringIO("")
    >>> io_repr(sio)
    'StringIOwith {name=None, mode=None, closed=False, readable=True, writable=True, is_textio=True, encoding=None, errors=None, line_buffering=False, newlines=None}'
    >>> sio.close()
    """
    if not is_io(fd):
        return repr(fd)

    repr_content = [
        f"name={repr(getattr(fd, 'name', None))}",
        f"mode={repr(getattr(fd, 'mode', None))}",
        f"closed={repr(fd.closed)}",
        f"readable={repr(fd.readable())}",
        f"writable={repr(fd.writable())}",
        f"is_textio={repr(is_textio(fd))}",
    ]
    if isinstance(fd, io.StringIO) or isinstance(fd, io.TextIOWrapper):
        repr_content.extend([
            f"encoding={repr(fd.encoding)}",
            f"errors={repr(fd.errors)}",
            f"line_buffering={repr(fd.line_buffering)}",
            f"newlines={repr(fd.newlines)}",
        ])

    if isinstance(fd, io.TextIOWrapper):
        repr_content.append(f"write_through={repr(fd.write_through)}")

    buffer: Optional[FDType] = None
    for possible_name in (
            "buffer",  # used by io.TextIOWrapper
            "_buffer",  # Used by lzma.LZMAFile
            "fileobj",  # Used by gzip.GzipFile
            "raw",  # Used by gzip._GZipReader
            "_fd",  # Used by IOProxy
            "_fp"  # Used by bz2.BZ2File
    ):
        try:
            buffer = getattr(fd, possible_name)
        except AttributeError:
            pass

    rets = fd.__class__.__name__ + "with {" + ", ".join(repr_content) + "}"
    if buffer is not None:
        rets += ", with buffer=<" + io_repr(buffer) + ">"
    return


def is_textio(fd: FDType) -> bool:
    """
    Determine whether a given IO is TextIO.
    A TextIO is usually an instance of :py:class:`io.TextIOBase` (:py:class:`io.TextIOWrapper`).
    Its :py:func:`read` function gnerates :py:class:`str`.
    A BinaryIO generates :py:class:`bytes` instead

    Example using :py:class:`io.StringIO`:

    >>> sio = io.StringIO("")
    >>> is_textio(sio)
    True
    >>> sio.close()

    Example using :py:class:`io.BytesIO`:

    >>> bio = io.BytesIO(b"")
    >>> is_textio(bio)
    False
    >>> bio.close()

    Example using File IO in text mode:

    >>> fio = open(__file__,"r")
    >>> is_textio(fio)
    True
    >>> fio.close()

    >>> fsio = open(__file__,"rt")
    >>> is_textio(fsio)
    True
    >>> fsio.close()

    Example using File IO in binary mode:

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

    Supports: :py:class:`typing.IO`, :py:class:`io.IOBase` and subclasses.

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

    >>> lzio = lzma.open(io.BytesIO(b"AAAA"))
    >>> is_io(lzio)
    True
    >>> lzio.close()

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
    Check whether a given object is a valid path.

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


@overload
def type_check(obj: FDType) -> Literal[True]:
    ...


@overload
def type_check(obj: PathType) -> Literal[False]:
    ...


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
    Determine line endings of a file descriptor. If failed, will return OS default.

    This accepts both binary and text IO while return type would always in string.

    :param fd: Input file descriptor.
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

    def __repr__(self):
        return f"IOProxy of IO <" + repr(self._fd) + ">"

    @property
    def mode(self) -> str:
        """
        Give the modestring of underlying file descriptor.

        .. warning ::
            This function not always give a string! See following situation:

            >>> import io, gzip
            >>> sio = io.StringIO()
            >>> sio.mode
            'NA'
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

    def write(self, s: bytes) -> int:  # type: ignore
        self._ensure_open()
        return self._fd.write(s)

    def writelines(self, lines: Iterable[bytes]) -> None:  # type: ignore
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


class TqdmReaderProxy(AbstractReader[AnyStr], Generic[AnyStr]):
    _tqdm: tqdm

    def __init__(self, fd: FDType):
        super().__init__(fd)
        self._tqdm = tqdm(
            desc=f"Reading from {self.name}",
            total=wc_c_io(self._fd),
            unit='B',
            unit_scale=True,
            unit_divisor=1024
        )

    def read(self, size: int = -1) -> AnyStr:
        self._ensure_open()
        update_bytes = self._fd.read(size)
        self._tqdm.update(len(update_bytes))
        return update_bytes

    def readline(self, limit: int = -1) -> AnyStr:
        self._ensure_open()
        update_bytes = self._fd.readline(limit)
        self._tqdm.update(len(update_bytes))
        return update_bytes

    def readlines(self, hint: int = -1) -> List[AnyStr]:
        self._ensure_open()
        update_bytes_arr = self._fd.readlines(hint)
        self._tqdm.update(sum(map(len, update_bytes_arr)))
        return update_bytes_arr


class LineReaderProxy(AbstractReader[AnyStr], Generic[AnyStr]):

    def read(self, *args, **kwargs):
        self._ensure_open()
        raise TypeError(f"Type {type(self._fd)} is Line Reader!")

    def readline(self, *args, **kwargs) -> AnyStr:
        self._ensure_open()
        return self._fd.readline(-1)  # Size limit canceled.

    def __iter__(self) -> Iterator[AnyStr]:
        """Iterator over lines with removal of trailing line feed (``LF``, ``\n``)/carriage return (``CR``, ``\r``)"""
        self._ensure_open()
        line = self.readline()
        if not line:
            return
        else:
            if isinstance(line, bytes):
                strip_content = b'\r\n'
            else:
                strip_content = '\r\n'
            yield line.rstrip(strip_content)
        while True:
            line = self.readline()
            if not line:
                break
            yield line.rstrip(strip_content)


class TqdmLineReaderProxy(LineReaderProxy[AnyStr], Generic[AnyStr]):
    _tqdm: tqdm

    def __init__(self, fd: FDType):
        super().__init__(fd)
        self._tqdm = tqdm(
            desc=f"Reading from {self.name}",
            total=wc_l_io(self._fd) + 1,
            unit='L',
            unit_scale=True,
            unit_divisor=1000
        )

    def readline(self, *args, **kwargs) -> AnyStr:
        self._ensure_open()
        update_bytes = self._fd.readline(-1)
        self._tqdm.update(1)
        return update_bytes


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


def file_open(  # type: ignore
        file_path: str,
        mode: ModeEnum = ModeEnum.READ,
        is_binary: bool = False,
        encoding: str = "UTF-8",
        compression: Optional[str] = "inferred",
        newline: Optional[str] = None,
        compression_level: int = 0,
        tqdm_reader_by_line: Optional[bool] = None,
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
            rfd = BinaryReaderProxy(rfd)
        else:
            rfd = TextReaderProxy(io.TextIOWrapper(
                rfd,
                encoding=encoding,
                newline=newline
            ))
        if tqdm_reader_by_line is True:
            return TqdmLineReaderProxy(rfd)
        elif tqdm_reader_by_line is False:
            return TqdmReaderProxy(rfd)
        else:
            return rfd

    elif mode == ModeEnum.WRITE:
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        open(file_path, mode="a").close()
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
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        open(file_path, mode="a").close()
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
) -> IOProxy[str]:
    ...


@overload
def get_reader(
        path_or_fd: PathOrFDType,
        is_binary: Literal[True],
        **kwargs
) -> IOProxy[bytes]:
    ...


def get_reader(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False, True] = False,
        **kwargs
) -> IOProxy:
    if type_check(path_or_fd):
        return wrap_io(path_or_fd)  # type: ignore
    else:
        path_or_fd = convert_path_to_str(path_or_fd)  # type: ignore
        return file_open(
            path_or_fd,
            mode=ModeEnum.READ,
            is_binary=is_binary,
            **kwargs
        )


@overload
def get_writer(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False],
        **kwargs
) -> IOProxy[str]:
    ...


@overload
def get_writer(
        path_or_fd: PathOrFDType,
        is_binary: Literal[True],
        **kwargs
) -> IOProxy[bytes]:
    ...


def get_writer(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False, True] = False,
        **kwargs
) -> IOProxy:
    if type_check(path_or_fd):
        return wrap_io(path_or_fd)  # type: ignore
    else:
        path_or_fd = convert_path_to_str(path_or_fd)  # type: ignore
        return file_open(
            path_or_fd,
            mode=ModeEnum.WRITE,
            is_binary=is_binary,
            **kwargs
        )


@overload
def get_appender(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False],
        **kwargs
) -> IOProxy[str]:
    ...


@overload
def get_appender(
        path_or_fd: PathOrFDType,
        is_binary: Literal[True],
        **kwargs
) -> IOProxy[bytes]:
    ...


def get_appender(
        path_or_fd: PathOrFDType,
        is_binary: Literal[False, True] = False,
        **kwargs
) -> IOProxy:
    if type_check(path_or_fd):
        return wrap_io(path_or_fd)  # type: ignore
    else:
        path_or_fd = convert_path_to_str(path_or_fd)  # type: ignore
        return file_open(
            path_or_fd,
            mode=ModeEnum.APPEND,
            is_binary=is_binary,
            **kwargs
        )
