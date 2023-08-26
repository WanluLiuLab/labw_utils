"""
``labw_utils.stdlib_helper.shutil_helper`` -- Enhanced :py:mod:`shutils` Module

More shell-like utilities.

.. versionadded:: 1.0.2
"""

__all__ = ("rm_rf", "wc_c", "wc_l", "wc_l_io", "wc_c_io", "touch")

import os
import shutil
import time

from labw_utils.commonutils.lwio import FDType, get_reader, wc_c_io, wc_l_io
from labw_utils.commonutils.lwio.file_system import get_abspath, file_exists, is_soft_link
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import IO, Callable, Optional

_lh = get_logger(__name__)


def readlink_f(path: str) -> str:
    """
    Remove soft links out of the path and return its absolute form,
    just like what is done by GNU CoreUtils ``readlink -f``.

    .. warning::
        Not tested, have previously-reported bugs.

    :param path: Input relative path
    :return: What you get from ``readlink -f``.

    .. versionadded:: 1.0.2
    """
    path = path.rstrip(os.sep)
    if path == "":
        return path
    while is_soft_link(path):
        path = get_abspath(os.path.realpath(path))
    return path


def wc_l(filename: str, opener: Optional[Callable[[str], IO]] = None) -> int:
    """
    Count lines in a file.

    :param filename: Input filename
    :param opener: Function to open this file. I.e., return an IO object.
    :return: Line number.

    .. versionadded:: 1.0.2
    """
    fd: IO
    if opener is None:
        fd = get_reader(filename, is_binary=True)
    else:
        fd = opener(filename)
    return wc_l_io(fd)


def wc_c(filename: str, opener: Optional[Callable[[str], IO]] = None) -> int:
    """
    Count the number of chars inside a file, i.e. File length.

    :param filename: Input filename
    :param opener: Function to open this file. I.e., return an IO object.
    :return: File length.

    .. versionadded:: 1.0.2
    """
    fd: FDType
    if opener is None:
        fd = get_reader(filename, is_binary=True)
    else:
        fd = opener(filename)
    return wc_c_io(fd)


def touch(
    filename: str, time_ns: Optional[float] = None, change_a_time: bool = True, change_m_time: bool = False
) -> None:
    """
    touch: ensure the existence of a file, just like GNU CoreUtils touch.

    .. seealso :: :manpage:`touch(1)`

    :param filename: The filename you wish to touch.
    :param change_m_time: Whether to change modification time.
    :param change_a_time: Whether to change access time.
    :param time_ns: Change time to this time instead of current time. Use :py:obj:`None` to use current time.
    :raises IsADirectoryError: If targeted filename is a directory.

    .. versionadded:: 1.0.2
    """
    filename = get_abspath(filename)
    if os.path.isdir(filename):
        raise IsADirectoryError(f"File '{filename}' is a directory")
    elif not file_exists(filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        open(filename, mode="a").close()
    if time_ns is None:
        time_ns = time.time_ns()
    if change_a_time:
        atime = time_ns
    else:
        atime = os.path.getatime(filename)
    if change_m_time:
        mtime = time_ns
    else:
        mtime = os.path.getmtime(filename)
    os.utime(filename, ns=(int(atime), int(mtime)))


def rm_rf(path: str) -> None:
    """
    Remove path recursively from the filesystem, just like rm -rf
    Should not complain on non-existing files.

    :param path: The path you wish to remove

    Compare with :py:func:`shutil.rmtree`,
    this function can remove files and things that do not exist without complaining.

    .. versionadded:: 1.0.2
    """
    dbg_head = "rm(path='" + path + "')"
    try:
        if os.path.isdir(path) and not os.path.islink(path):
            _lh.debug(f"{dbg_head} is a directory")
            shutil.rmtree(path)
        elif os.path.exists(path):
            _lh.debug(f"{dbg_head} is a file")
            os.remove(path)
    except FileNotFoundError:
        _lh.debug(f"{dbg_head} not exist")
