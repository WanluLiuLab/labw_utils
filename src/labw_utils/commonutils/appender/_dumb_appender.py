from labw_utils.commonutils.appender import BaseTableAppender, TableAppenderConfig
from labw_utils.typing_importer import Any, Tuple, Sequence


class DumbTableAppender(BaseTableAppender):
    """Appender that does nothing."""

    def __init__(self, filename: str, header: Sequence[str], tac: TableAppenderConfig):
        super().__init__(filename, header, tac)

    def _get_real_filename_hook(self):
        return ""

    def _create_file_hook(self):
        """Not needed"""
        pass

    def append(self, body: Sequence[Any]):
        """Not needed"""
        pass

    def close(self):
        """Not needed"""
        pass

    def flush(self):
        """Not needed"""
        pass
