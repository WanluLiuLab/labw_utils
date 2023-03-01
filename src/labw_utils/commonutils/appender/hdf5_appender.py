import os


from labw_utils import UnmetDependenciesError

try:
    import tables as pt
except ImportError:
    raise UnmetDependenciesError("pytables")

try:
    import pandas as pd
except ImportError:
    raise UnmetDependenciesError("pandas")

from labw_utils.commonutils.appender._pandas_table_appender import PandasDictBufferAppender


class HDF5TableAppender(PandasDictBufferAppender):
    def _get_n_lines_actually_written_hook(self) -> int:
        df = pd.read_hdf(self._real_filename, key="df")
        if isinstance(df, pd.DataFrame):
            return df.shape[0]
        else:
            raise TypeError

    def _get_real_filename_hook(self):
        self._real_filename = ".".join((self.filename, "hdf5"))

    def _create_file_hook(self):
        """Function not required"""
        pass

    def _write_hook(self, df: pd.DataFrame):
        df.to_hdf(self._real_filename, key="df", format='table', append=os.path.exists(self._real_filename))
