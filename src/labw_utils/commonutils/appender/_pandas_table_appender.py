from abc import ABC

from labw_utils import UnmetDependenciesError

try:
    import pandas as pd
except ImportError:
    raise UnmetDependenciesError("pandas")


from labw_utils.commonutils.appender.typing import BaseDictBufferAppender


class PandasDictBufferAppender(BaseDictBufferAppender, ABC):

    def flush(self) -> pd.DataFrame:
        df = pd.DataFrame.from_dict(data=self._buff)
        return df
