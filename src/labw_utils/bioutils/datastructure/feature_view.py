"""
gtf_view -- Ordered list of GTF records.
"""
from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.feature import Feature
from labw_utils.commonutils.lwio.file_system import should_regenerate
from labw_utils.commonutils.stdlib_helper import pickle_helper
from labw_utils.typing_importer import Iterable, Iterator, List, Sequence, SequenceProxy

FVPKL_VERSION = "1.0.0"


class FeatureView(Iterable[Feature]):
    _l: List[Feature]

    def __init__(
            self,
            l: Iterable[Feature],
            shortcut: bool = False
    ):
        if not shortcut:
            self._l = list(l)
        else:
            self._l = l  # type: ignore

    @property
    def l(self) -> Sequence[Feature]:
        return SequenceProxy(self._l)

    @classmethod
    def from_gtf(
            cls,
            gtf_file_path: str,
            show_tqdm: bool = True
    ):
        index_file_path = f"{gtf_file_path}.{FVPKL_VERSION}.fvpkl.xz"
        if not should_regenerate(gtf_file_path, index_file_path):
            try:
                return cls.from_fvpkl(index_file_path)
            except TypeError:
                pass

        new_instance = cls(
            list(GtfIterator(gtf_file_path, show_tqdm=show_tqdm)),
            shortcut=True
        )
        new_instance.to_fvpkl(index_file_path)
        return new_instance

    @classmethod
    def from_fvpkl(
            cls,
            index_file_path: str
    ):
        (fvpkl_version, new_instance) = pickle_helper.load(index_file_path)
        if fvpkl_version == FVPKL_VERSION:
            return new_instance
        else:
            raise TypeError(f"Version {fvpkl_version} (file) and {FVPKL_VERSION} (library) mismatch")

    def to_fvpkl(self, index_file_path: str):
        pickle_helper.dump((FVPKL_VERSION, self), index_file_path)

    def __iter__(self) -> Iterator[Feature]:
        yield from self._l
