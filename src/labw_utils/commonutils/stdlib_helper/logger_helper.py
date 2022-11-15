"""
logger_helper.py -- System-Wide Logger.

Features
--------

It can also read environment variable named 'LOG_LEVEL' and fallback to DEBUG (10) by default.

Usage
-----

There are also :py:func:`set_level` and :py:func:`get_logger`,
which is only snake-case wrappers for those contents inside :py:mod:`logging` standard module.
"""

import logging
import os
from logging import DEBUG, WARNING, ERROR, FATAL
from typing import Union

__all__ = (
    'DEBUG',
    'WARNING',
    'ERROR',
    'FATAL',
    'set_level',
    'get_logger'
)

_lh = logging.getLogger()


def _get_formatter(level: int) -> logging.Formatter:
    if level > logging.DEBUG:
        log_format = '%(asctime)s\t[%(levelname)s] %(message)s'
    else:
        log_format = '%(asctime)s %(name)s:%(lineno)d::%(funcName)s\t[%(levelname)s]\t%(message)s'
    return logging.Formatter(log_format)


def set_level(level: Union[str, int], quiet: bool = True) -> int:
    """
    Set the global logging level, and update the format.
    The log will be more verbose if the level is below debug.

    # FIXME: File set DEBUG.
    :raise ValueError: Raise this error if level not exist
    """
    _lh.setLevel(level)
    this_level = _lh.getEffectiveLevel()
    if not quiet:
        _lh.info(f'Resetting log level: {logging.getLevelName(this_level)}')

    logging.basicConfig(
        # level=this_level,
        handlers=[
            logging.StreamHandler()
        ]
    )
    file_handler = logging.FileHandler(filename="log.log")
    file_handler.setLevel(logging.DEBUG)
    # logging.root.setLevel(logging.DEBUG)
    for handler in logging.root.handlers:
        handler.formatter = _get_formatter(this_level)
    logging.root.addHandler(file_handler)
    return this_level


if "_global_level" not in locals() and "_global_level" not in globals():
    # set the global log-level.
    # Will read from LOG_LEVEL environment variable.
    # And fall to DEBUG if fails.
    _global_level = os.environ.get('LOG_LEVEL')
    if _global_level is None:
        _global_level = logging.INFO
    _global_level = set_level(_global_level)


def get_logger(name: str):
    """
    A Simple logging.getLogger() wrapper.
    """
    return logging.getLogger(name)
