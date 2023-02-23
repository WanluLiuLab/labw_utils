import importlib
from typing import Type, Iterator, Tuple

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils.appender.typing import BaseTableAppender

POSSIBLE_APPENDER_PATHS = (
    "labw_utils.commonutils.appender.tsv_appender",
    "labw_utils.commonutils.appender.lzmatsv_appender",
    "labw_utils.commonutils.appender.lz77tsv_appender",
    "labw_utils.commonutils.appender.dumb_appender",
    "labw_utils.commonutils.appender.hdf5_appender",
    "labw_utils.commonutils.appender.parquet_appender",
    "labw_utils.commonutils.appender.sqlite3_appender",
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


def load_table_appender_class(name: str) -> Type[BaseTableAppender]:
    """
    Return a known tracer.

    :raise UnmetDependenciesError: If dependencies are unmet
    """
    for possible_path in POSSIBLE_APPENDER_PATHS:
        try:
            mod = importlib.import_module(possible_path)
            return getattr(mod, name)
        except (ModuleNotFoundError, AttributeError):
            continue
    raise ModuleNotFoundError


def list_table_appender() -> Iterator[Tuple[str, str]]:
    for possible_path in POSSIBLE_APPENDER_PATHS:
        try:
            mod = importlib.import_module(possible_path)

            for k, v in mod.__dict__.items():
                if k.__contains__("Appender") and not k.__contains__("Base"):
                    try:
                        yield k, v.__doc__.strip().splitlines()[0]
                    except AttributeError:
                        yield k, "No docs available"
        except (ModuleNotFoundError, UnmetDependenciesError):
            continue
