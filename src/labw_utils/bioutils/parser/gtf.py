"""
TODO: docs

.. versionadded:: 1.0.2
"""
from labw_utils.bioutils.parser import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.record.feature import FeatureInterface, DEFAULT_GTF_QUOTE_OPTIONS
from labw_utils.bioutils.record.gtf import parse_record, format_string, GTFParsingError
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.typing_importer import Iterable, Iterator


class GtfIterator(BaseFileIterator, Iterable[FeatureInterface]):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    .. versionchanged:: 1.0.3
        `onerror` was added to skip erroneous lines when parsing malformed GTF
    """

    filetype: str = "GTF"
    record_type = FeatureInterface
    onerror: str

    def __init__(self, filename: str, onerror: str = "alert", **kwargs):
        """
        :param filename: Filename to be parsed.
        :param onerror: What to do on error. Use ``skip`` to skip the erroneous line.
        """
        super().__init__(filename, **kwargs)
        self.onerror = onerror

    def __iter__(self) -> Iterator[FeatureInterface]:
        for line in get_tqdm_line_reader(self.filename):
            if line.startswith("#") or line == "":
                continue
            try:
                yield parse_record(line)
            except GTFParsingError as e:
                if self.onerror == "skip":
                    pass
                else:
                    raise e


class GtfIteratorWriter(BaseIteratorWriter):
    """
    TODO: docs

    .. versionadded:: 1.0.2
    """

    filetype: str = "GTF"
    record_type = FeatureInterface

    def __init__(self, filename: str, quote: str = DEFAULT_GTF_QUOTE_OPTIONS, **kwargs):
        super().__init__(filename, **kwargs)
        self._fd = get_writer(self._filename)
        self._quote = quote

    @staticmethod
    def write_iterator(
        iterable: Iterable[FeatureInterface],
        filename: str,
        prefix_annotations: Iterable[str] = None,
        quote: str = DEFAULT_GTF_QUOTE_OPTIONS,
        **kwargs,
    ):
        with GtfIteratorWriter(filename, quote) as writer:
            if prefix_annotations is not None:
                for annotation in prefix_annotations:
                    writer.write_comment(annotation)
            for feature in iterable:
                writer.write(feature)

    def write(self, record: FeatureInterface) -> None:
        self._fd.write(format_string(record, quote=self._quote) + "\n")

    def write_comment(self, comment: str):
        self._fd.write("#" + comment + "\n")
