import multiprocessing
import os
from abc import abstractmethod, ABC
from multiprocessing import synchronize
from typing import List, Any, Dict


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
    _header: List[str]
    _real_filename: str
    _tac: TableAppenderConfig

    @property
    def filename(self) -> str:
        return self._filename

    @property
    def header(self) -> List[str]:
        return self._header

    @property
    def real_filename(self) -> str:
        return self._real_filename

    def __init__(self, filename: str, header: List[str], tac: TableAppenderConfig):
        self._filename = filename
        self._header = header
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
    def append(self, body: List[Any]) -> None:
        raise NotImplementedError

    @abstractmethod
    def close(self) -> None:
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class BaseDictBufferAppender(BaseTableAppender, ABC):
    _h0: str
    _buff: Dict[str, List[Any]]
    _write_mutex: synchronize.Lock
    _buff_mutex: synchronize.Lock

    def __init__(self, filename: str, header: List[str], tac: TableAppenderConfig):
        super().__init__(filename, header, tac)
        self._buff_mutex = multiprocessing.Lock()
        self._write_mutex = multiprocessing.Lock()
        self._buff = {}
        self._h0 = self.header[0]

    def append(self, body: List[Any]):
        with self._buff_mutex:
            if self._buff == {}:
                self._buff = dict(zip(self.header, map(lambda x: [x], body)))
            else:
                for header_item, body_item in zip(self.header, body):
                    self._buff[header_item].append(body_item)
            if len(self) == self._tac.buffer_size:
                df = self._flush()
                self._buff = {}
                with self._write_mutex:
                    self._write_hook(df)

    @abstractmethod
    def _write_hook(self, df: Any):
        raise NotImplementedError

    @abstractmethod
    def _flush(self) -> Any:
        raise NotImplementedError

    def __len__(self):
        if self._buff == {}:
            return 0
        return len(self._buff[self._h0])

    def flush(self):
        if len(self) == 0:
            return
        df = self._flush()
        self._buff = {}
        with self._write_mutex:
            self._write_hook(df)

    def close(self):
        self.flush()
