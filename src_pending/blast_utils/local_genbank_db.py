import glob
import io
import os.path
import re
import shutil

from Bio import Entrez

from labw_utils.commonutils.lwio import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Iterable, Dict, Tuple

_lh = get_logger(__name__)

_GB_VERSION_REGEX = re.compile(r"^VERSION\s+(\S+)$")


def parse_gb_simple(s: io.TextIOWrapper) -> Iterable[Tuple[str, str]]:
    """VERSION and full contents"""
    current_record_version = ""
    current_full_str = ""
    for l in s:
        l = l.strip("\n\r")
        current_full_str += l
        current_full_str += "\n"
        if l == "//":
            yield current_record_version, current_full_str
            current_record_version = ""
            current_full_str = ""
            continue
        version_match = _GB_VERSION_REGEX.match(l)
        if version_match is not None:
            current_record_version = version_match.group(1)


class LocalGenBankDB:
    _path: str
    _cache: Dict[str, str]

    @property
    def path(self):
        return self._path

    def __init__(self, path: str):
        self._path = os.path.abspath(path)
        os.makedirs(self._path, exist_ok=True)
        self._cache = {}
        for p in glob.glob(os.path.join(self._path, "*.gbk")):
            dbid = os.path.basename(p).replace(".gbk", "")
            self._cache[dbid] = p

    def prepare(self, gbaccs: Iterable[str]):
        gbaccs = list(filter(lambda s: s not in self._cache, set(gbaccs)))
        _lh.info("Preparation of %d accessions required", len(gbaccs))
        if len(gbaccs) == 0:
            return
        with io.StringIO() as tmpf:
            handle = Entrez.efetch(
                "protein",
                id=",".join(gbaccs),
                rettype="gb",
                retmode="text",
            )
            shutil.copyfileobj(handle, tmpf)
            handle.close()
            tmpf.seek(0)

            for version, contents in parse_gb_simple(tmpf):
                gb_path = os.path.join(self.path, version + ".gbk")
                with get_writer(gb_path, is_binary=False) as w:
                    w.write(contents)
                self._cache[version] = gb_path

    def get(self, gbacc: str) -> str:
        return self._cache[gbacc]
