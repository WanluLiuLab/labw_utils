from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Mapping, Any, Optional, Callable, Dict
from labw_utils import UnmetDependenciesError

from labw_utils.stdlib.cpy311 import tomllib

try:
    import tomli_w
except ImportError:
    raise UnmetDependenciesError("tomli_w")


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
        path: str,
        title: str,
        validate_versions: Optional[Callable[[Mapping[str, Any]], None]] = None
) -> Dict[str, Any]:
    """
    Read and validate TOML files with metadata.
    """
    with open(path, "rb") as reader:
        in_dict = tomllib.load(reader)
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
