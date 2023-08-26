"""
``labw_utils.commonutils.appender`` -- Appenders of relational data format

.. warning::
    is not finished -- do not use.
"""

import importlib
import os
from abc import ABC, abstractmethod
import pkgutil

from labw_utils import UnmetDependenciesError
from labw_utils.stdlib.cpy310.pkgutil import resolve_name
from labw_utils.typing_importer import Type, Iterator, Tuple, Any, Sequence

_POSSIBLE_APPENDER_PATHS = tuple(
    f"{__package__}.{spec.name}"
    for spec in pkgutil.iter_modules(resolve_name(__package__).__spec__.submodule_search_locations)
)


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
    """
    .. versionadded:: 1.0.2

    .. versionchanged:: 1.0.3
        List instances changed from :py:class:`tuple` to :py:class:`list`.
    """

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

    def __init__(self, filename: str, header: Sequence[str], tac: TableAppenderConfig):
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
    def append(self, body: Sequence[Any]) -> None:
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
    for possible_path in _POSSIBLE_APPENDER_PATHS:
        try:
            mod = importlib.import_module(possible_path)
            return getattr(mod, name)
        except (ModuleNotFoundError, AttributeError, UnmetDependenciesError, IndexError):
            continue
    raise ModuleNotFoundError


def list_table_appender() -> Iterator[Tuple[str, str]]:
    """
    List table appenders that can be imported and their documentations.
    """
    models = []
    for possible_path in _POSSIBLE_APPENDER_PATHS:
        try:
            mod = importlib.import_module(possible_path)

            for k, v in mod.__dict__.items():
                if (
                    k.__contains__("Appender")
                    and not k.__contains__("Base")
                    and not k.__contains__("Config")
                    and k not in models
                ):
                    try:
                        yield k, v.__doc__.strip().splitlines()[0]
                    except AttributeError:
                        yield k, "No docs available"
                    models.append(k)
        except (ModuleNotFoundError, AttributeError, UnmetDependenciesError, ImportError):
            continue
