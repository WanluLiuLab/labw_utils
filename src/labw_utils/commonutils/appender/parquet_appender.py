import os

from labw_utils import UnmetDependenciesError

try:
    import fastparquet as fp
except ImportError:
    raise UnmetDependenciesError("fastparquet")

try:
    import pandas as pd
except ImportError:
    raise UnmetDependenciesError("pandas")

from labw_utils.commonutils.appender._pandas_table_appender import PandasDictBufferAppender


class ParquetTableAppender(PandasDictBufferAppender):

    def _get_n_lines_actually_written_hook(self) -> int:
        return pd.read_parquet(self._real_filename).shape[0]

    def _get_real_filename_hook(self):
        self._real_filename = ".".join((self.filename, "parquet"))

    def _create_file_hook(self):
        """Function not needed"""
        pass

    def _write_hook(self, df: pd.DataFrame):
        fp.write(self._real_filename, df, append=os.path.exists(self._real_filename))
