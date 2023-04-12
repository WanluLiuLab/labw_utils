"""
rule_based_ioproxy -- An IO implementation which may read/write archives or texts.

This IO wrapper will open archives of ``.gz``, ``.xz``, ``.lzma`` or ``.bz2`` suffix,
with standard IO functions like ``read``,
with iteration and context management (aka., ``__enter__`` or ``with`` statement) support.
"""

from __future__ import annotations

__all__ = (
    "determine_line_endings",
    "get_reader",
    "get_appender",
    "get_writer",
    "RuleBasedSequentialReader"
)

import labw_utils.commonutils.io as cio
from labw_utils.commonutils.io import PathOrFDType, IOProxy
from labw_utils.commonutils.io import determine_line_endings as raw_determine_line_endings
from labw_utils.commonutils.io.rules import FileRuleRing

_ILLEGAL_READONLY_IO_EXCEPTION = OSError("Illegal operation on Sequential Read-Only IO")


def determine_line_endings(file_path_or_fd: PathOrFDType) -> str:
    """
    Determine file endings using :py:class:`RuleBasedIOProxy`.
    """
    return raw_determine_line_endings(RuleBasedIOProxy(file_path_or_fd))


class RuleBasedIOProxy(IOProxy):

    def __init__(
            self,
            path_or_fd: PathOrFDType,
            *args,
            **kwargs
    ):
        """
        :param path_or_fd: The filename to be opened.
          If is file descriptor, will be passed to IOProxy.
        :param args: Positional arguments passed to underlying opener.
        :param kwargs: Keyword arguments passed to underlying opener.
        """
        if cio.type_check(path_or_fd):
            super().__init__(path_or_fd)
        else:
            super().__init__(FileRuleRing.open(path_or_fd, *args, **kwargs))


class RuleBasedSequentialReader(RuleBasedIOProxy):
    """
    A sequential reader that is based on :py:class:`RuleBasedIOProxy`
    which does not support random-access functions like :py:func:`seek`
    or writing functions like :py:func:`write`
    """

    def seek(self, *args, **kwargs) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise _ILLEGAL_READONLY_IO_EXCEPTION

    def seekable(self) -> bool:
        """False"""
        return False

    def truncate(self, *args, **kwargs) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise _ILLEGAL_READONLY_IO_EXCEPTION

    def writable(self) -> bool:
        """False"""
        return False

    def write(self, *args, **kwargs) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise _ILLEGAL_READONLY_IO_EXCEPTION

    def writelines(self, *args, **kwargs) -> None:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise _ILLEGAL_READONLY_IO_EXCEPTION


def get_reader(path_or_fd: PathOrFDType, is_binary: bool = False, **kwargs) -> IOProxy:
    """
    Rule-based :py:func:`labw_utils.commonutils.io.get_reader`.

    Will fall through to :py:class:`IOProxy` if input is file descriptor.
    """
    if cio.type_check(path_or_fd):
        return IOProxy(path_or_fd)
    else:
        if is_binary:
            mode = "rb"
        else:
            mode = "rt"
        return RuleBasedIOProxy(path_or_fd, mode=mode, **kwargs)


def get_writer(path_or_fd: PathOrFDType, is_binary: bool = False, **kwargs) -> IOProxy:
    """
    Rule-based :py:func:`labw_utils.commonutils.io.get_writer`.

    Will fall through to :py:class:`IOProxy` if input is file descriptor.
    """
    if cio.type_check(path_or_fd):
        return IOProxy(path_or_fd)
    else:
        if is_binary:
            mode = "wb"
        else:
            mode = "wt"
        return RuleBasedIOProxy(path_or_fd, mode=mode, **kwargs)


def get_appender(path_or_fd: PathOrFDType, is_binary: bool = False, **kwargs) -> IOProxy:
    """
    Rule-based :py:func:`labw_utils.commonutils.io.get_appender`.

    Will fall through to :py:class:`IOProxy` if input is file descriptor.
    """
    if cio.type_check(path_or_fd):
        return IOProxy(path_or_fd)
    else:
        if is_binary:
            mode = "ab"
        else:
            mode = "at"
        return RuleBasedIOProxy(path_or_fd, mode=mode, **kwargs)
