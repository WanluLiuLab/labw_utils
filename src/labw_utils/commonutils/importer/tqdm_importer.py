"""
tqdm_importer.py -- Import `tqdm` without messing up stderr

This module imports `tqdm`, the progress bar implementation on Python.

If import is failed or stderr is not a Pseudo Terminal,
will use a home-made fallback which is more silent.
"""
import os
import sys

__all__ = (
    "tqdm"
)

from labw_utils.commonutils.importer import _silent_tqdm

if os.getenv("SPHINX_BUILD") is not None:
    try:
        import tqdm as _external_tqdm
    except ImportError:
        _external_tqdm = None
else:
    _external_tqdm = None

if sys.stderr.isatty() and _external_tqdm is not None:
    tqdm = _external_tqdm.tqdm
else:
    tqdm = _silent_tqdm.tqdm
