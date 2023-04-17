"""
safe_io -- A Safe Wrapper for :py:mod:`commonutils.io`

This is a "safe" IO. It does follow things:

On reader, it ensures existence of file being read by throwing errors.

On writer or appender, it ensures existence of file by :py:func:`touch`-ing it first.
If the argument is an IO, will also check whether it is writable.
"""

__all__ = (
    "get_reader",
    "get_writer",
    "get_appender"
)

import labw_utils.commonutils.io as cio
import labw_utils.commonutils.io.rule_based_ioproxy as rio
from labw_utils.commonutils.io import PathOrFDType
from labw_utils.commonutils.io import file_system, IOProxy
from labw_utils.commonutils.stdlib_helper import shutil_helper


def get_reader(path_or_fd: PathOrFDType, **kwargs) -> IOProxy:
    """
    Safe rule-based :py:func:`labw_utils.commonutils.io.get_reader`.

    Will fall through to :py:class:`IOProxy` if input is file descriptor.
    """
    if cio.type_check(path_or_fd):
        return IOProxy(path_or_fd)
    if file_system.file_exists(path_or_fd, allow_special_paths=True):
        return rio.get_reader(path_or_fd, **kwargs)
    else:
        raise FileNotFoundError(f"File {path_or_fd} not found!")


def get_writer(path_or_fd: PathOrFDType, **kwargs) -> IOProxy:
    """
    Safe rule-based :py:func:`labw_utils.commonutils.io.get_writer`.

    Will fall through to :py:class:`IOProxy` if input is file descriptor.
    """
    if cio.type_check(path_or_fd):
        if path_or_fd.writable():
            raise TypeError("Attempt to write on read-only IO")
        return IOProxy(path_or_fd)
    if not file_system.file_exists(path_or_fd, allow_special_paths=True):
        shutil_helper.touch(path_or_fd)
    return rio.get_writer(path_or_fd, **kwargs)


def get_appender(path_or_fd: PathOrFDType, **kwargs) -> IOProxy:
    """
    Safe rule-based :py:func:`labw_utils.commonutils.io.get_appender`.

    Will fall through to :py:class:`IOProxy` if input is file descriptor.
    """
    if cio.type_check(path_or_fd):
        if path_or_fd.writable():
            raise TypeError("Attempt to write on read-only IO")
        return IOProxy(path_or_fd)
    if file_system.file_exists(path_or_fd, allow_special_paths=True):
        return rio.get_appender(path_or_fd, **kwargs)
    else:
        shutil_helper.touch(path_or_fd)
