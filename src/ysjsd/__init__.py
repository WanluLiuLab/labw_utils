import platform
import sys
import time
from typing import Mapping, Any, Final

import flask

import labw_utils

__version__ = "0.1.0"

from labw_utils.mlutils.io_helper import AbstractTOMLSerializable


class YSJSDConfig(AbstractTOMLSerializable):
    _title: Final[str] = "ysjsd"
    _name: str
    _description: str
    _ysjs_port: str
    _var_directory: str
    _config_file_path: str

    @staticmethod
    def _validate_versions(versions: Mapping[str, Any]) -> None:
        pass

    @staticmethod
    def _dump_versions() -> Mapping[str, Any]:
        return {
                "labw_utils": labw_utils.__version__,
                "ysjsd": __version__,
                "flask": flask.__version__,
                "python": ".".join(map(str, (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)))
            }

    @staticmethod
    def _dump_metadata() -> Mapping[str, Any]:
        return {
            "time_gmt": time.asctime(time.gmtime()),
            "time_local": time.asctime(time.localtime()),
            "platform": platform.uname()._asdict()
        }

    def __init__(
            self,
            name: str,
            description: str,
            ysjs_port: str,
            var_directory: str,
            config_file_path: str
    ):
        self._name = name
        self._description = description
        self._ysjs_port = ysjs_port
        self._var_directory = var_directory
        self._config_file_path = config_file_path

    @property
    def ysjs_port(self) -> str:
        return self._ysjs_port

    def to_dict(self) -> Mapping[str, str]:
        return {
            "name": self._name,
            "description": self._description,
            "ysjs_port": self._ysjs_port,
            "var_directory": self._var_directory,
            "config_file_path": self._config_file_path
        }

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

