import multiprocessing
from abc import ABC, abstractmethod
from multiprocessing import synchronize

from labw_utils.commonutils.appender import BaseTableAppender, TableAppenderConfig
from labw_utils.typing_importer import Dict, Any, Tuple, List, Sequence


class BaseDictBufferAppender(BaseTableAppender, ABC):
    _h0: str
    """Name of the first column."""
    _buff: Dict[str, List[Any]]
    """Column-oriented buffer."""
    _write_mutex: synchronize.Lock
    _buff_mutex: synchronize.Lock

    def __init__(self, filename: str, header: Sequence[str], tac: TableAppenderConfig):
        super().__init__(filename, header, tac)
        self._buff_mutex = multiprocessing.Lock()
        self._write_mutex = multiprocessing.Lock()
        self._init_buffer()
        self._h0 = self.header[0]

    def _init_buffer(self):
        self._buff = {k: [] for k in self.header}

    def append(self, body: Sequence[Any]):
        with self._buff_mutex:
            for header_item, body_item in zip(self.header, body):
                self._buff[header_item].append(body_item)
            if len(self) == self._tac.buffer_size:
                self.flush()

    @abstractmethod
    def _write_hook(self, df: Any):
        raise NotImplementedError

    @abstractmethod
    def _convert_dict_to_df(self) -> Any:
        """
        Convert the buffer dict to intermediate format.
        """
        raise NotImplementedError

    def __len__(self):
        return len(self._buff[self._h0])

    def flush(self):
        if len(self) == 0:
            return
        with self._write_mutex:
            self._write_hook(self._convert_dict_to_df())
        self._init_buffer()

    def close(self):
        self.flush()
