from abc import ABC

from labw_utils import UnmetDependenciesError

try:
    import pandas as pd
except ImportError:
    raise UnmetDependenciesError("pandas")

from labw_utils.commonutils.appender._base_dict_buffer_appender import BaseDictBufferAppender


class BasePandasDictBufferAppender(BaseDictBufferAppender, ABC):
    def _convert_dict_to_df(self) -> pd.DataFrame:
        df = pd.DataFrame.from_dict(data=self._buff)
        return df
