import functools
import glob
import shutil
import subprocess
from labw_utils.typing_importer import Optional

import pandas as pd

from blast_utils.ncbi_taxdb import NON_EXIST, TaxonDB, CELLULAR_ORGANISMS
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Sequence

_lh = get_logger(__name__)


class AbstractConfig:
    def which(self, what: str) -> Optional[str]:
        return shutil.which(what)

    def require(self, program_names: Sequence[str]):
        for program_name in program_names:
            if self.which(program_name) is None:
                raise RuntimeError(f"Program {program_name} is missing!")

def merge_table(src_path_glob: str, dst_path: str):
    _lh.info("Merging table files...")
    functools.reduce(
        lambda a, b: pd.concat((a, b)),
        map(lambda p: pd.read_csv(p, sep="\t"), glob.glob(src_path_glob)),
    ).to_csv(
        path_or_buf=dst_path,
        sep="\t",
        index=False,
    )
