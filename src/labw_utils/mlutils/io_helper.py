"""
``labw_utils.mlutils.io_helper`` -- Compressed serialization formats.

Here provides compressed readers and writers for Numpy and serialization formats,
which can significantly reduce disk size.

We also have an abstract base class that allows programmers to create their own configuration class.

The compression algorithm would be Lempel-Ziv Markov Chain Algorithm (LZMA) version 2 used in
`7-Zip <https://www.7-zip.org>`_. The implementation is provided Python standard library :py:mod:`lzma`.

.. warning::
    Since Python's standard LZMA implementation is single-threaded,
    it might be extremely slow when compressing large objects!

.. versionadded:: 1.0.0

.. versionchanged:: 1.0.1
    Serialization of Torch :py:class:`torch.Tensor` was removed as it is compressed natively.
"""

from __future__ import annotations

__all__ = ("read_np_xz", "write_np_xz")

import lzma

import numpy as np
import numpy.lib.format as npy_format
import numpy.typing as npt


def read_np_xz(path: str) -> npt.NDArray:
    """
    Reader of compressed Numpy serialization format

    .. versionadded:: 1.0.0
    """
    with lzma.open(path, "rb") as reader:
        return npy_format.read_array(reader)


def write_np_xz(array: npt.NDArray, path: str) -> None:
    """
    Writer of compressed Numpy serialization format

    .. versionadded:: 1.0.0
    """
    with lzma.open(path, "wb", preset=9) as writer:
        npy_format.write_array(writer, np.asanyarray(array))
