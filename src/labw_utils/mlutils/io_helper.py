"""
Here provides compressed readers and writers for Pickle, Numpy and Torch serialization formats,
which can significantly reduce disk size.

We also have an abstract base class that allows programmers to create their own configuration class.

The compression algorithm would be Lempel-Ziv Markov Chain Algorithm (LZMA) version 2 used in
`7-Zip <https://www.7-zip.org>`_. The implementation is provided Python standard library :py:mod:`lzma`.

.. warning::
    Since Python's standard LZMA implementation is single-threaded,
    it might be extremely slow to compress large objects!
"""

from __future__ import annotations

__all__ = (
    "read_np_xz",
    "read_tensor_xz",
    "write_np_xz",
    "write_tensor_xz"
)

import lzma
from typing import Any, Union, Mapping

import numpy as np
import numpy.lib.format as npy_format
import numpy.typing as npt

try:
    import torch
except ImportError:
    torch = None


def read_np_xz(path: str) -> npt.NDArray:
    """Reader of compressed Numpy serialization format"""
    with lzma.open(path, "rb") as reader:
        return npy_format.read_array(reader)


def write_np_xz(array: npt.NDArray, path: str) -> None:
    """Writer of compressed Numpy serialization format"""
    with lzma.open(path, "wb", preset=9) as writer:
        npy_format.write_array(writer, np.asanyarray(array))


if torch is not None:
    def read_tensor_xz(path: str) -> Union[torch.Tensor, Mapping[str, Any], torch.nn.Module]:
        """Reader of compressed Torch serialization format"""
        with lzma.open(path, "rb") as reader:
            return torch.load(reader)


    def write_tensor_xz(array: Union[torch.Tensor, Mapping[str, Any], torch.nn.Module], path: str) -> None:
        """Writer of compressed Torch serialization format"""
        with lzma.open(path, "wb", preset=9) as writer:
            torch.save(array, writer)
