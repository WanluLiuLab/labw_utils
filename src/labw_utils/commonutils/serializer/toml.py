"""
``labw_utils.commonutils.serializer.toml`` -- TOML serialization interfaces.

Is not finished -- do not use.

>>> from labw_utils.typing_importer import Final
>>> import io
>>> class A(AbstractTOMLSerializable):
...     _a: int
...     _b: int
...     _title: Final[str] = "A"
...
...     def __init__(self, a: int, b: int):
...         self._a = a
...         self._b = b
...
...     def to_dict(self) -> Mapping[str, Any]:
...         return {"a": self._a, "b": self._b}
...
...     @classmethod
...     def from_dict(cls, in_dict: Mapping[str, Any]):
...         return cls(**in_dict)
...
...     @staticmethod
...     def _dump_versions() -> Optional[Mapping[str, Any]]:
...         return {"version": 1}
...
...     @staticmethod
...     def _dump_metadata() -> Optional[Mapping[str, Any]]:
...         return {}
...
...     @staticmethod
...     def _validate_versions(versions: Mapping[str, Any]) -> None:
...         return None
>>> sio = io.BytesIO()
>>> A(1, 2).save(sio)
>>> print(str(sio.getvalue(), encoding="UTF-8"))
[A]
a = 1
b = 2
<BLANKLINE>
[version_info]
version = 1
<BLANKLINE>
[metadata]
<BLANKLINE>

.. versionadded:: 1.0.2
"""

from __future__ import annotations

from abc import ABC, abstractmethod
import os

# TODO: Finish this module
if os.getenv("SPHINX_BUILD") is not None:
    __all__ = []
else:
    __all__ = [
        "TOMLRepresentableInterface",
        "AbstractTOMLSerializable",
        "read_toml_with_metadata",
        "write_toml_with_metadata",
    ]

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils import lwio
from labw_utils.stdlib.cpy311 import tomllib
from labw_utils.typing_importer import Any, Optional, Callable
from labw_utils.typing_importer import Mapping

try:
    import pytest

    tomli_w = pytest.importorskip("tomli_w")

except ImportError:
    try:
        import tomli_w
    except ImportError as e:
        raise UnmetDependenciesError("tomli_w") from e

from labw_utils.commonutils.serializer import SerializableInterface


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
    path_of_fd: lwio.PathOrFDType,
    title: str,
    validate_versions: Optional[Callable[[Mapping[str, Any]], None]] = None,
) -> Mapping[str, Any]:
    """
    Read and validate TOML files with metadata.
    """
    if (
        lwio.is_io(path_of_fd)
        and path_of_fd.readable()
        and not lwio.is_textio(path_of_fd)
    ):
        in_dict = tomllib.load(path_of_fd)
    elif lwio.is_path(path_of_fd):
        with lwio.file_open(
            file_path=path_of_fd, mode=lwio.ModeEnum.READ, is_binary=True
        ) as reader:
            in_dict = tomllib.load(reader)
    else:
        raise TypeError(f"Type {type(path_of_fd)} not supported!")
    if "version_info" in in_dict and validate_versions is not None:
        validate_versions(in_dict.pop("version_info"))
    if "metadata" in in_dict:
        _ = in_dict.pop("metadata")
    return in_dict.pop(title)


def write_toml_with_metadata(
    obj: Mapping[str, Any],
    title: str,
    path_of_fd: lwio.PathOrFDType,
    dump_versions: Optional[Callable[[], Optional[Mapping[str, Any]]]] = None,
    dump_metadata: Optional[Callable[[], Optional[Mapping[str, Any]]]] = None,
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

    if (
        lwio.is_io(path_of_fd)
        and path_of_fd.writable()
        and not lwio.is_textio(path_of_fd)
    ):
        tomli_w.dump(retd, path_of_fd)
    elif lwio.is_path(path_of_fd):
        with lwio.file_open(
            file_path=path_of_fd, mode=lwio.ModeEnum.WRITE, is_binary=True
        ) as writer:
            tomli_w.dump(retd, writer)
    else:
        raise TypeError(f"Type {type(path_of_fd)} not supported!")


class AbstractTOMLSerializable(SerializableInterface, TOMLRepresentableInterface, ABC):
    """
    Abstract Base Class of something that can be represented as TOML.

    Should be used as configuration class.
    """

    _title: str

    @classmethod
    def load(cls, path: str, **kwargs):
        return cls.from_dict(
            read_toml_with_metadata(path, cls._title, cls._validate_versions)
        )

    def save(self, path_of_fd: lwio.PathOrFDType, **kwargs) -> None:
        write_toml_with_metadata(
            self.to_dict(),
            self._title,
            path_of_fd,
            self._dump_versions,
            self._dump_metadata,
        )
