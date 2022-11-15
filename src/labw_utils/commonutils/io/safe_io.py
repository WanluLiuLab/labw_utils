"""
safe_io -- A Safe Wrapper for :py:mod:`commonutils.io`
======================================================

This is a "safe" IO. It does follow things:

On reader, it ensures existence of file being read by throwing errors.

On writer or appender, it ensures existence of file by :py:func:`touch`-ing it first.
"""

from labw_utils.commonutils.io import file_system
from labw_utils.commonutils.io import get_appender as _get_appender
from labw_utils.commonutils.io import get_reader as _get_reader
from labw_utils.commonutils.io import get_writer as _get_writer
from labw_utils.commonutils.io._ioproxy import IOProxyType
from labw_utils.commonutils.stdlib_helper import shutil_helper
from labw_utils.devutils.decorators import copy_doc

__all__ = (
    "get_reader",
    "get_writer",
    "get_appender"
)

from labw_utils.commonutils.typing import PathOrFDType


@copy_doc(_get_reader)
def get_reader(path_or_fd: PathOrFDType, **kwargs) -> IOProxyType:
    if file_system.file_exists(path_or_fd, allow_special_paths=True):
        return _get_reader(path_or_fd, **kwargs)
    else:
        raise FileNotFoundError(f"File {path_or_fd} not found!")


@copy_doc(_get_writer)
def get_writer(path_or_fd: PathOrFDType, **kwargs) -> IOProxyType:
    if not file_system.file_exists(path_or_fd, allow_special_paths=True):
        shutil_helper.touch(path_or_fd)
    return _get_writer(path_or_fd, **kwargs)


@copy_doc(_get_appender)
def get_appender(path_or_fd: PathOrFDType, **kwargs) -> IOProxyType:
    if file_system.file_exists(path_or_fd, allow_special_paths=True):
        return _get_appender(path_or_fd, **kwargs)
    else:
        shutil_helper.touch(path_or_fd)
