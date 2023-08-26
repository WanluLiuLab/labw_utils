"""
tqdm_reader.py -- Reader with Progress Bar

Here are wrappings for basic IO classes & functions in
:py:mod:`labw_utils.commonutils.lwio` with additional progress bar.
"""

__all__ = ("get_tqdm_reader", "get_tqdm_line_reader")

import labw_utils.commonutils.lwio as cio


def get_tqdm_reader(*args, **kwargs) -> cio.IOProxy:
    return cio.get_reader(*args, **kwargs, tqdm_reader=True)


def get_tqdm_line_reader(*args, **kwargs) -> cio.IOProxy:
    return cio.get_reader(*args, **kwargs, tqdm_reader=True, line_reader=True)
