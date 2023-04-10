"""
labw_utils.devutils.sphinx_helper -- Helpers of Sphinx documentation system
"""
__all__ = (
    "convert_ipynb_to_myst",
)

import glob
import os

import jupytext

from labw_utils.commonutils.io.file_system import should_regenerate
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def convert_ipynb_to_myst(
        source_dir: str
):
    for fn in glob.glob(os.path.join(source_dir, "**", "*.ipynb"), recursive=True):
        fn = os.path.abspath(fn)
        if "ipynb_checkpoints" in fn or "_build" in fn:
            continue
        dst_fn = fn + ".md"
        if should_regenerate(fn, dst_fn):
            _lh.info(f"CONVERT {fn} -> {dst_fn}: START")
            try:
                jupytext.write(
                    nb=jupytext.read(fn, fmt="notebook"),
                    fp=dst_fn,
                    fmt="md:myst"
                )
            except Exception:
                _lh.warning(f"CONVERT {fn} -> {dst_fn}: ERR")
            _lh.info(f"CONVERT {fn} -> {dst_fn}: FIN")
        else:
            _lh.info(f"CONVERT {fn} -> {dst_fn}: REFUSE TO OVERWRITE NEWER FILE")

