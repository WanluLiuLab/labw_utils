import gzip

from labw_utils.commonutils.appender._tsv_appender import TSVTableAppender


class LZ77TSVTableAppender(TSVTableAppender):
    """Append to GZipped TSV."""

    def _get_real_filename_hook(self) -> str:
        return ".".join((self.filename, "tsv", "gz"))

    def _create_file_hook(self):
        with gzip.open(self._real_filename, mode="wt") as writer:
            writer.write("\t".join(self.header) + "\n")

    def _write_hook(self, df: str):
        with gzip.open(self._real_filename, mode="at") as writer:
            writer.write(df)
