from typing import List, Any

from labw_utils.commonutils.appender import BaseTableAppender, TableAppenderConfig


class DumbTableAppender(BaseTableAppender):
    def __init__(self, filename: str, header: List[str], tac: TableAppenderConfig):
        super().__init__(filename, header, tac)

    def _get_real_filename_hook(self):
        return ""

    def _create_file_hook(self):
        """Not needed"""
        pass

    def append(self, body: List[Any]):
        """Not needed"""
        pass

    def close(self):
        """Not needed"""
        pass

    def flush(self):
        """Not needed"""
        pass
