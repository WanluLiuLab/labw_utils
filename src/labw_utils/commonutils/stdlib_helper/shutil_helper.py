"""
shell_utils -- Enhanced :py:mod:`shutils` Module

More shell-like utilities.
"""

import gzip
import os
import shutil
from typing import IO, Callable

from labw_utils.commonutils.io import get_reader
from labw_utils.commonutils.io.file_system import get_abspath, file_exists, is_soft_link


def readlink_f(path: str) -> str:
    """
    Remove soft links out of the path and return its abstract form, just like what is done by GNU CoreUtils readlink -f.

    Can be used to trace symlink of symlink.

    :param path: Input relative path
    :return: What you get from readlink -f
    """
    path = path.rstrip(os.sep)
    if path == '':
        return path
    while is_soft_link(path):
        path = get_abspath(os.path.realpath(path))
    return path


def wc_l(filename: str, opener: Callable[[str], IO] = None) -> int:
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


def wc_c(filename: str, opener: Callable[[str], IO] = None) -> int:
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


def wc_l_io(fd: IO) -> int:
    """
    Count lines in a file.

    :param fd: An IO object.
    :return: Line number.
    """
    if fd.seekable():
        curr_pos = fd.tell()
    else:
        return -1
    reti = 0
    fd.seek(0)
    segment = fd.read(1)
    if len(segment) == 0:
        return 0
    else:
        if isinstance(segment, str):
            reti += segment.count("\n")
            while True:
                segment = fd.read(1024)
                if len(segment) == 0:
                    break
                reti += segment.count("\n")
        else:
            while True:
                reti += segment.count("\n")
                segment = fd.read(1024)
                if len(segment) == 0:
                    break
                reti += segment.count(b"\n")
    if fd.tell() != 0 and reti == 0:
        reti = 1  # To keep similar behaviour to GNU WC
    fd.seek(curr_pos)
    return reti


def wc_c_io(fd: IO) -> int:
    """
    Count the number of chars inside a file, i.e. File length.

    :param fd: An IO object.
    :return: File length.
    """
    if fd.seekable():
        curr_pos = fd.tell()
        fd.seek(0, 2)
        reti = fd.tell()
        fd.seek(curr_pos)
        return reti
    else:
        return -1


def touch(filename: str) -> None:
    """
    touch: ensure the existence of a file, just like GNU CoreUtils touch.
    Please note that this version of touch CANNOT change file attributes like times.

    :param filename: The filename you wish to touch
    """
    filename = get_abspath(filename)
    if file_exists(filename):
        return
    elif os.path.isdir(filename):
        raise IsADirectoryError(f"File '{filename}' is a directory")
    else:
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        open(filename, mode="a").close()


def rm_rf(path: str) -> None:
    """
    Remove path recursively from the filesystem, just like rm -rf
    Should not complain on non-existing files.

    :param path: The path you wish to remove

    Compare with :py:func:`shutil.rmtree`, this function can remove files.
    """
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def gz_compress(in_file: str, out_file: str, keep_in_file: bool = False, compresslevel: int = 9) -> None:
    with open(in_file, 'rb') as f_in, gzip.open(out_file, 'wb', compresslevel=compresslevel) as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_in_file:
        rm_rf(in_file)


def gz_decompress(in_file: str, out_file: str, keep_in_file: bool = False) -> None:
    with gzip.open(in_file, 'rb') as f_in, open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_in_file:
        rm_rf(in_file)
