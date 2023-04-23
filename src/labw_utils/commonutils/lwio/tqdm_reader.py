"""
tqdm_reader.py -- Reader with Progress Bar

Here are wrappings for basic IO classes & functions in
:py:mod:`labw_utils.commonutils.lwio` with additional progress bar.
"""

__all__ = (
    "TqdmReader",
    "TqdmLineReader",
    "get_tqdm_reader",
    "get_tqdm_line_reader"
)

import labw_utils.commonutils.lwio as cio
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio import get_reader, AbstractReader
from labw_utils.commonutils.stdlib_helper import shutil_helper
from labw_utils.devutils.decorators import copy_doc
from labw_utils.typing_importer import Iterator, AnyStr, List


class _BaseTqdmReader(AbstractReader):
    _tqdm: tqdm

    def __enter__(self, *args, **kwargs):
        self._tqdm.__enter__(*args, **kwargs)
        return self

    def __exit__(self, *args, **kwargs):
        return self._tqdm.__exit__(*args, **kwargs)


class TqdmReader(_BaseTqdmReader):
    """
    A very simple tqdm reader with :py:func:``read`` and :py:func:``readline`` functions.
    """

    @copy_doc(RuleBasedIOProxy.__init__)
    def __init__(self, filename: str, *args, **kwargs):
        super().__init__(filename, *args, **kwargs)
        self._tqdm = tqdm(
            desc=f"Reading {filename}", total=shutil_helper.wc_c_io(self._fd), unit='B', unit_scale=True,
            unit_divisor=1024
        )

    @copy_doc(RuleBasedIOProxy.readline)
    def readline(self, *args, **kwargs) -> AnyStr:
        update_bytes = super().readline(*args, **kwargs)
        self._tqdm.update(len(update_bytes))
        return update_bytes

    @copy_doc(RuleBasedIOProxy.read)
    def read(self, size: int = -1) -> AnyStr:
        update_bytes = super().read(size)
        self._tqdm.update(len(update_bytes))
        return update_bytes

    @copy_doc(RuleBasedIOProxy.readlines)
    def readlines(self, *args, **kwargs) -> List[AnyStr]:
        update_bytes_arr = super().readlines(*args, **kwargs)
        self._tqdm.update(sum(map(len, update_bytes_arr)))
        return update_bytes_arr


class TqdmLineReader(_BaseTqdmReader):
    """
    A very simple tqdm reader with only :py:func:`readline` functions.

    .. warning::
        This class has different :py:func:`__iter__` method!
        It would remove trailing line feeds.
    """

    @copy_doc(RuleBasedIOProxy.__init__)
    def __init__(self, filename: str, *args, **kwargs):
        super().__init__(filename, *args, **kwargs)
        self._tqdm = tqdm(
            desc=f"Reading {filename}",
            total=shutil_helper.wc_l_io(self._fd) + 1,
            unit='L'
        )

    def read(self, *args, **kwargs):
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError('Illegal operation read.')

    def readlines(self, *args, **kwargs):
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError('Illegal operation read.')

    @copy_doc(RuleBasedIOProxy.readline)
    def readline(self, *args, **kwargs) -> AnyStr:
        update_bytes = super().readline(-1)  # Size limit canceled.
        self._tqdm.update(1)
        return update_bytes

    def __iter__(self) -> Iterator[str]:
        """Iterator over lines with removal of trailing line feed (``LF``, ``\n``)/carriage return (``CR``, ``\r``)"""
        while True:
            line = self.readline()
            if not line:
                break
            yield line.rstrip('\n\r')


@copy_doc(get_reader)
def get_tqdm_reader(path_or_fd: cio.PathOrFDType, is_binary: bool = False, **kwargs) -> RuleBasedSequentialReader:
    """
    :py:func:`get_reader`-like wrapper for :py:class:`TqdmReader`.

    If input is IO, will pass through to :py:class`RuleBasedSequentialReader`.
    """
    if cio.type_check(path_or_fd):
        return RuleBasedSequentialReader(path_or_fd)
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return TqdmReader(path_or_fd, mode=mode, **kwargs)


def get_tqdm_line_reader(path_or_fd: cio.PathOrFDType, is_binary: bool = False, **kwargs) -> RuleBasedSequentialReader:
    """
    :py:func:`get_reader`-like wrapper for :py:class:`TqdmLineReader`.

    If input is IO, will pass through to :py:class`RuleBasedSequentialReader`.
    """
    if cio.type_check(path_or_fd):
        return RuleBasedSequentialReader(path_or_fd)
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return TqdmLineReader(path_or_fd, mode=mode, **kwargs)
