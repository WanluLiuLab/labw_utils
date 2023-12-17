"""
``labw_utils.commonutils.importer.tqdm_importer`` -- Import ``tqdm`` without messing up stderr

This module imports `tqdm <https://tqdm.github.io/>`_, the progress bar implementation in Python.

If import is failed or stderr is not a Pseudo Terminal,
will use a home-made fallback which is more silent.

Supports :envvar:`TQDM_IMPL` environment variable to override the defaults.

List of environment variables:

.. envvar:: TQDM_IMPL

Allowed values:

- ``EXTERNAL``: To use official ``tqdm``.
- ``SILENT``: Use internal one.
- others: Auto decide.

.. versionadded:: 1.0.2
"""
import os
import sys

__all__ = ("tqdm",)

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
    # Force use internal when building docs
    _external_tqdm = None

if os.getenv("TQDM_IMPL") == "EXTERNAL":
    pass
elif os.getenv("TQDM_IMPL") == "SILENT":
    _external_tqdm = None
else:
    if not sys.stderr.isatty():
        _external_tqdm = None

if _external_tqdm is not None:
    tqdm = _external_tqdm.tqdm
else:
    tqdm = _silent_tqdm.tqdm
