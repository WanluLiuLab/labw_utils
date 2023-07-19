import lzma

from labw_utils.commonutils.appender._tsv_appender import TSVTableAppender


class LZMATSVTableAppender(TSVTableAppender):
    """Append to XZipped TSV."""

    def _get_real_filename_hook(self) -> str:
        return ".".join((self.filename, "tsv", "xz"))

    def _create_file_hook(self):
        with lzma.open(self._real_filename, mode="wt") as writer:
            writer.write("\t".join(self.header) + "\n")

    def _write_hook(self, df: str):
        with lzma.open(self._real_filename, mode="at") as writer:
            writer.write(df)
