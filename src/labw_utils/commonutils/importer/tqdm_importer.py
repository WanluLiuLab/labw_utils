"""
tqdm_importer.py -- Import `tqdm` without messing up stderr

This module imports `tqdm`, the progress bar implementation on Python.

If import is failed or stderr is not a Pseudo Terminal,
will use a home-made fallback which is more silent.
"""
import os
import sys

__all__ = (
    "tqdm",
)

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.importer import _silent_tqdm

_lh = get_logger(__name__)

if os.getenv("SPHINX_BUILD") is None:
    try:
        import tqdm as _external_tqdm
    except ImportError:
        _external_tqdm = None
        _lh.warning("Import official tqdm failed! will use builtin instead")
else:
    _external_tqdm = None
    _lh.warning("Sphinx environment detected!")

if os.getenv("TQDM_IMPL") == "EXTERNAL":
    pass
elif os.getenv("TQDM_IMPL") == "SILENT":
    _external_tqdm = None
else:
    if not sys.stderr.isatty():
        _lh.warning("STDERR is not TTY!")
        _external_tqdm = None

if _external_tqdm is not None:
    tqdm = _external_tqdm.tqdm
else:
    tqdm = _silent_tqdm.tqdm
