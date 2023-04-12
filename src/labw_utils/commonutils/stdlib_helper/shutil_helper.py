"""
labw_utils.stdlib_helper.shutil_helper -- Enhanced :py:mod:`shutils` Module

More shell-like utilities.
"""

__all__ = (
    "rm_rf",
    "wc_c",
    "wc_l",
    "wc_l_io",
    "wc_c_io",
    "touch"
)

import os
import shutil
import time
from typing import IO, Callable, Optional

from labw_utils.commonutils.io.rule_based_ioproxy import get_reader
from labw_utils.commonutils.io.file_system import get_abspath, file_exists, is_soft_link
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def readlink_f(path: str) -> str:
    """
    Remove soft links out of the path and return its absolute form,
    just like what is done by GNU CoreUtils ``readlink -f``.

    Can be used to trace symlink of symlink.

    TODO: to be re-checked.

    :param path: Input relative path
    :return: What you get from ``readlink -f``.
    """
    path = path.rstrip(os.sep)
    if path == '':
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
    """
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
    """
    if opener is None:
        fd = get_reader(filename, is_binary=True)
    else:
        fd = opener(filename)
    return wc_c_io(fd)


def wc_l_io(fd: IO, block_size: int = 4096) -> int:
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
                reti += block.count("\n")
        else:
            reti += block.count(b"\n")
            while True:
                block = fd.read(block_size)
                if len(block) == 0:
                    break
                reti += block.count(b"\n")
    if fd.tell() != 0 and reti == 0:
        reti = 1  # To keep similar behaviour to GNU WC
    fd.seek(curr_pos)
    return reti


def wc_c_io(fd: IO) -> int:
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


def touch(
        filename: str,
        time_ns: Optional[float] = None,
        change_a_time: bool = True,
        change_m_time: bool = False
) -> None:
    """
    touch: ensure the existence of a file, just like GNU CoreUtils touch.

    :param filename: The filename you wish to touch.
    :param change_m_time: Whether to change modification time.
    :param change_a_time: Whether to change access time.
    :param time_ns: Change time to this time instead of current time. Use :py:obj:`None` to use current time.
    :raises IsADirectoryError: If targeted filename is a directory.
    """
    filename = get_abspath(filename)
    if file_exists(filename):
        pass
    elif os.path.isdir(filename):
        raise IsADirectoryError(f"File '{filename}' is a directory")
    else:
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

    Compare with :py:func:`shutil.rmtree`, this function can remove files.
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
