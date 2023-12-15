"""
file_system.py -- Basic Filesystem Functions

Here are very low-level filesystem functions used by other Python modules,
like :py:mod:`commonutils.lwio.safe_io` or :py:mod:`commonutils.stdlib_helper.shutil_helper`.
"""

__all__ = (
    "get_abspath",
    "file_exists",
    "directory_exists",
    "is_soft_link",
    "should_regenerate",
    "to_safe_filename",
)

import os
import stat

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def to_safe_filename(filename: str) -> str:
    """
    Get safe filename on any platform
    """

    filename = filename[:256]
    filename = (
        filename.replace(" ", "_")
        .replace("\t", "_")
        .replace("\\", "_")
        .replace(":", "_")
        .replace("*", "_")
        .replace("?", "_")
        .replace('"', "_")
        .replace("<", "_")
        .replace(">", "_")
        .replace("|", "_")
        .replace("/", "_")
    )
    if filename.endswith(".") or filename.endswith(" "):
        filename += "_"
    if filename.upper() in (
        "CON",
        "PRN",
        "AUX",
        "CLOCK$",
        "NUL",
        "COM1",
        "COM2",
        "COM3",
        "COM4",
        "COM5",
        "COM6",
        "COM7",
        "COM8",
        "COM9",
        "LPT1",
        "LPT2",
        "LPT3",
        "LPT4",
        "LPT5",
        "LPT6",
        "LPT7",
        "LPT8",
        "LPT9",
    ):
        filename += "_"
    return filename


def get_abspath(path: str) -> str:
    """
    Get the absolute path of a given relative path. Can be a non-existent file or directory.

    The default one provided by os.path.abspath() or os.path.realpath()
    cannot perform tide expand (That is, expand ~ to your HOME directory).

    :param path: The relative path
    :return: The absolute path
    """
    if path == "":
        return path
    path = os.path.abspath(os.path.expanduser(path))
    return path


def file_exists(path: str, allow_special_paths: bool = True) -> bool:
    """
    Check the existence of a file. Can check regular and special files.

    :param path: The path you wish to examine.
    :param allow_special_paths: Whether this is a special file, may be block device, etc.
    """
    if allow_special_paths:
        return os.path.exists(path) and not os.path.isdir(path)
    else:
        return os.path.isfile(path)


def directory_exists(path: str) -> bool:
    """
    Check the existence of a directory. Can check regular and special files.

    :param path: The path you wish to examine.
    """
    return os.path.exists(path) and os.path.isdir(path)


def is_soft_link(path: str) -> bool:
    """
    This function checks whether provided path is a soft (symbolic) link.

    :param path: The path you wish to examine.
    """
    return stat.S_ISLNK(os.stat(path, follow_symlinks=False)[0])


def should_regenerate(src_path: str, dst_path: str, verbose: bool = True) -> bool:
    """
    This function assumes ``dst_path`` is generated from ``src_path``
    and determine whether re-generation is necessary.

    :param src_path: The source path.
    :param dst_path: The destination path.
    :param verbose: Whether to be verbose. Will print logs at INFO level.
    """
    if file_exists(dst_path):
        if os.path.getmtime(dst_path) - os.path.getmtime(src_path) > 0:
            return False
        elif verbose:
            _lh.info("REGEN %s -> %s: New older than source", src_path, dst_path)
            return True
    elif verbose:
        _lh.info("REGEN %s -> %s: File does not exist", src_path, dst_path)
        return True
