"""
safe_io -- A Safe Wrapper for :py:mod:`commonutils.lwio`

.. deprecated:: 1.0.1
    Use :py:func:`labw_utils.commonutils.lwio.file_open` instead.
"""

__all__ = ("get_reader", "get_writer", "get_appender")

from labw_utils.commonutils import lwio


def get_reader(*args, **kwargs) -> lwio.IOProxy:
    """
    Deprecated. Do not use.
    """
    return lwio.get_reader(*args, **kwargs)


def get_writer(*args, **kwargs) -> lwio.IOProxy:
    """
    Deprecated. Do not use.
    """
    return lwio.get_writer(*args, **kwargs)


def get_appender(*args, **kwargs) -> lwio.IOProxy:
    """
    Deprecated. Do not use.
    """
    return lwio.get_appender(*args, **kwargs)
