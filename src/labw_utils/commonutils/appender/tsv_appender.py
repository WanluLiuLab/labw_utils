from labw_utils.commonutils.appender.typing import BaseDictBufferAppender
from labw_utils.commonutils.shell_utils import wc_l


class TSVTableAppender(BaseDictBufferAppender):

    def _get_real_filename_hook(self):
        self._real_filename = ".".join((self.filename, "tsv"))

    def flush(self) -> str:
        return "\n".join(map(
            lambda x: "\t".join(map(repr, x)),  # x is [COLUMN]
            zip(*self._buff.values())
        )) + "\n"

    def _create_file_hook(self):
        with open(self._real_filename, mode="wt") as writer:
            writer.write("\t".join(self.header) + "\n")

    def _write_hook(self, df: str):
        with open(self._real_filename, mode="at") as writer:
            writer.write(df)

    def _get_n_lines_actually_written_hook(self) -> int:
        return wc_l(self._real_filename) - 1
