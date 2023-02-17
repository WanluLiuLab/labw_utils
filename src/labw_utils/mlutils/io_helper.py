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
    "write_tensor_xz",
    "SerializableInterface",
    "TOMLRepresentableInterface",
    "AbstractTOMLSerializable",
    "write_toml_with_metadata",
    "read_toml_with_metadata"
)

import lzma
import pickle
from abc import abstractmethod, ABC
from typing import Any, Union, Mapping, Dict, Callable, Optional

import numpy as np
import numpy.lib.format as npy_format
import numpy.typing as npt
import tomli
import tomli_w

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


    def write_tensor_xz(array: Union[torch.Tensor, Mapping[str, Any], nn.Module], path: str) -> None:
        """Writer of compressed Torch serialization format"""
        with lzma.open(path, "wb", preset=9) as writer:
            torch.save(array, writer)


class SerializableInterface(ABC):
    """
    Something that can be saved or loaded to files.
    """

    @classmethod
    def load(cls, path: str, **kwargs):
        """
        Load configuration from a file.

        :param path: Filename to read from.
        :return: New instance of corresponding class.
        """
        raise NotImplementedError

    def save(self, path: str, **kwargs) -> None:
        """
        Save the class contents with metadata.

        :param path: Filename to write to.
        """
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def _dump_versions() -> Optional[Mapping[str, Any]]:
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def _dump_metadata() -> Optional[Mapping[str, Any]]:
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def _validate_versions(versions: Mapping[str, Any]) -> None:
        raise NotImplementedError


class TOMLRepresentableInterface(ABC):
    """
    Interface of something that can be represented as TOML.

    Should be used as configuration class.
    """

    @abstractmethod
    def to_dict(self) -> Mapping[str, Any]:
        """Dump the item to a dictionary"""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        """Load the item from a dictionary"""
        raise NotImplementedError


def read_toml_with_metadata(
        path: str,
        title: str,
        validate_versions: Optional[Callable[[Mapping[str, Any]], None]] = None
) -> Dict[str, Any]:
    """
    Read and validate TOML files with metadata.
    """
    with open(path, "rb") as reader:
        in_dict = tomli.load(reader)
    if "version_info" in in_dict and validate_versions is not None:
        validate_versions(in_dict.pop("version_info"))
    if "metadata" in in_dict:
        _ = in_dict.pop("metadata")
    return in_dict.pop(title)


def write_toml_with_metadata(
        obj: Mapping[str, Any],
        title: str,
        path: str,
        dump_versions: Optional[Callable[[], Mapping[str, Any]]] = None,
        dump_metadata: Optional[Callable[[], Mapping[str, Any]]] = None
) -> None:
    """
    Write TOML files with metadata.
    """
    retd = {title: obj}
    if dump_versions is not None:
        version_info = dump_versions()
        if version_info is not None:
            retd["version_info"] = version_info
    if dump_metadata is not None:
        metadata = dump_metadata()
        if metadata is not None:
            retd["metadata"] = metadata
    with open(path, 'wb') as writer:
        tomli_w.dump(retd, writer)


class AbstractTOMLSerializable(
    SerializableInterface,
    TOMLRepresentableInterface,
    ABC
):
    """
    Abstract Base Class of something that can be represented as TOML.

    Should be used as configuration class.
    """
    _title: str

    @classmethod
    def load(cls, path: str, **kwargs):
        return cls.from_dict(read_toml_with_metadata(
            path,
            cls._title,
            cls._validate_versions
        ))

    def save(self, path: str, **kwargs) -> None:
        write_toml_with_metadata(
            self.to_dict(),
            self._title,
            path,
            self._dump_versions,
            self._dump_metadata
        )
