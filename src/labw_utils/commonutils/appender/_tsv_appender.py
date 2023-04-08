from labw_utils.commonutils.appender._base_dict_buffer_appender import BaseDictBufferAppender


class TSVTableAppender(BaseDictBufferAppender):

    def _get_real_filename_hook(self) -> str:
        return ".".join((self.filename, "tsv"))

    def _flush(self) -> str:
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