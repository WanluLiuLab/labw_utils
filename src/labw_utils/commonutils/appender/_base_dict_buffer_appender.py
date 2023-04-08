import multiprocessing
from abc import ABC, abstractmethod
from multiprocessing import synchronize
from typing import Dict, List, Any

from labw_utils.commonutils.appender import BaseTableAppender, TableAppenderConfig


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
