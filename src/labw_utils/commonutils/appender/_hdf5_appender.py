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

from labw_utils.commonutils.appender._base_pandas_dict_table_appender import BasePandasDictBufferAppender


class HDF5TableAppender(BasePandasDictBufferAppender):

    def _get_real_filename_hook(self) -> str:
        return ".".join((self.filename, "hdf5"))

    def _create_file_hook(self):
        """Function not required"""
        pass

    def _write_hook(self, df: pd.DataFrame):
        df.to_hdf(self._real_filename, key="df", format='table', append=os.path.exists(self._real_filename))
