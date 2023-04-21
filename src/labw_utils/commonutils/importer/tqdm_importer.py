"""
tqdm_importer.py -- Import `tqdm` without messing up stderr

This module imports `tqdm`, the progress bar implementation on Python.

If import is failed or stderr is not a Pseudo Terminal,
will use a home-made fallback which is more silent.
"""
import os
import sys

from labw_utils.typing_importer import Type

__all__ = (
    "tqdm",
)

tqdm: Type

if os.getenv("SPHINX_BUILD") is not None:
    from labw_utils.commonutils.importer._silent_tqdm import tqdm as tqdm
else:
    try:
        import tqdm as _external_tqdm
    except ImportError:
        _external_tqdm = None
    from labw_utils.commonutils.importer._silent_tqdm import tqdm as _silent_tqdm

    if sys.stderr.isatty() and _external_tqdm is not None:
        tqdm = _external_tqdm.tqdm
    else:
        tqdm = _silent_tqdm
