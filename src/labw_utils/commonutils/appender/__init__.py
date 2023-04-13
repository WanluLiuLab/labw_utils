"""
labw_utils.commonutils.appender -- Appenders of relational data format

TODO: Docs
"""

import importlib
import os
from abc import ABC, abstractmethod
from typing import Type, Iterator, Tuple, Any

from labw_utils import UnmetDependenciesError

POSSIBLE_APPENDER_PATHS = (
    "labw_utils.commonutils.appender._tsv_appender",
    "labw_utils.commonutils.appender._lzmatsv_appender",
    "labw_utils.commonutils.appender._lz77tsv_appender",
    "labw_utils.commonutils.appender._dumb_appender",
    "labw_utils.commonutils.appender._hdf5_appender",
    "labw_utils.commonutils.appender._parquet_appender",
    "labw_utils.commonutils.appender._sqlite3_appender",
)

AVAILABLE_TABLE_APPENDERS = {
    "DumbTableAppender": 'DumbTableAppender',
    "TSVTableAppender": 'TSVTableAppender',
    "LZMATSVTableAppender": 'LZMATSVTableAppender',
    "LZ77TSVTableAppender": 'LZ77TSVTableAppender',
    "HDF5TableAppender": 'HDF5TableAppender',
    "ParquetTableAppender": 'ParquetTableAppender',
    "SQLite3TableAppender": 'SQLite3TableAppender'
}


class TableAppenderConfig:
    _buffer_size: int
    """
    Buffering strategy. 1 for no buffering.
    """

    def __init__(self, buffer_size: int = 1):
        self._buffer_size = buffer_size

    @property
    def buffer_size(self) -> int:
        return self._buffer_size


class BaseTableAppender(ABC):
    _filename: str
    _header: Tuple[str, ...]
    _real_filename: str
    _tac: TableAppenderConfig

    @property
    def filename(self) -> str:
        return self._filename

    @property
    def header(self) -> Tuple[str, ...]:
        return self._header

    @property
    def real_filename(self) -> str:
        return self._real_filename

    def __init__(self, filename: str, header: Tuple[str, ...], tac: TableAppenderConfig):
        self._filename = filename
        self._header = tuple(header)
        self._real_filename = self._get_real_filename_hook()
        self._tac = tac
        if os.path.exists(self._real_filename):
            os.remove(self._real_filename)
        self._create_file_hook()

    @abstractmethod
    def _get_real_filename_hook(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def _create_file_hook(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def flush(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def append(self, body: Tuple[Any, ...]) -> None:
        raise NotImplementedError

    @abstractmethod
    def close(self) -> None:
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def load_table_appender_class(name: str) -> Type[BaseTableAppender]:
    """
    Return a known tracer.

    :raise UnmetDependenciesError: If dependencies are unmet.
    :raise ModuleNotFoundError: If is not found.
    """
    for possible_path in POSSIBLE_APPENDER_PATHS:
        try:
            mod = importlib.import_module(possible_path)
            return getattr(mod, name)
        except (ModuleNotFoundError, AttributeError):
            continue
    raise ModuleNotFoundError


def list_table_appender() -> Iterator[Tuple[str, str]]:
    """
    List table appenders that can be imported and their documentations.
    """
    models = []
    for possible_path in POSSIBLE_APPENDER_PATHS:
        try:
            mod = importlib.import_module(possible_path)

            for k, v in mod.__dict__.items():
                if k.__contains__("Appender") and not k.__contains__("Base") and not k.__contains__("Config") and not k in models:
                    try:
                        yield k, v.__doc__.strip().splitlines()[0]
                    except AttributeError:
                        yield k, "No docs available"
                    models.append(k)
        except (ModuleNotFoundError, UnmetDependenciesError):
            continue
