from labw_utils.commonutils.appender._base_dict_buffer_appender import BaseDictBufferAppender
from labw_utils.commonutils.lwio.safe_io import get_writer


class TSVTableAppender(BaseDictBufferAppender):
    """Append to TSV."""

    def _get_real_filename_hook(self) -> str:
        return ".".join((self.filename, "tsv"))

    def _convert_dict_to_df(self) -> str:
        return "\n".join(map(lambda x: "\t".join(map(repr, x)), zip(*self._buff.values()))) + "\n"  # x is [COLUMN]

    def _create_file_hook(self):
        with get_writer(self._real_filename, is_binary=False) as writer:
            writer.write("\t".join(self.header) + "\n")

    def _write_hook(self, df: str):
        with open(self._real_filename, mode="at") as writer:
            writer.write(df)
