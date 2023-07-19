import sqlite3

from labw_utils import UnmetDependenciesError

try:
    import pandas as pd
except ImportError:
    raise UnmetDependenciesError("pandas")

from labw_utils.commonutils.appender._base_pandas_dict_table_appender import BasePandasDictBufferAppender


class SQLite3TableAppender(BasePandasDictBufferAppender):
    """Append to SQLite. Requires Pandas."""

    def _get_real_filename_hook(self) -> str:
        return ".".join((self.filename, "sqlite3"))

    def _create_file_hook(self):
        """Not needed"""
        pass

    def _write_hook(self, df: pd.DataFrame):
        with sqlite3.connect(self._real_filename) as con:
            df.to_sql(name="db", con=con, if_exists="append", index=False)
