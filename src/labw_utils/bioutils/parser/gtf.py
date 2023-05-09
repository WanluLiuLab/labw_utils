from labw_utils.bioutils.parser import BaseFileIterator, BaseIteratorWriter
from labw_utils.bioutils.record.feature import FeatureInterface, DEFAULT_GTF_QUOTE_OPTIONS
from labw_utils.bioutils.record.gtf import parse_record, format_string
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.typing_importer import Iterable, Iterator


class GtfIterator(BaseFileIterator, Iterable[FeatureInterface]):
    filetype: str = "GTF"
    record_type = FeatureInterface

    def __iter__(self) -> Iterator[FeatureInterface]:
        for line in get_tqdm_line_reader(self.filename):
            if line.startswith('#') or line == '':
                continue
            yield parse_record(line)


class GtfIteratorWriter(BaseIteratorWriter):
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
            quote: str = DEFAULT_GTF_QUOTE_OPTIONS
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
        self._fd.write('#' + comment + "\n")
