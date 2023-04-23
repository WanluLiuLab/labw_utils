import json
import os
import shutil

from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

setup_basic_logger()

_lh = get_logger()

if __name__ == "__main__":
    with get_reader("fa_spec.json") as reader:
        fa_spec = json.load(reader)
    for fa_entry in fa_spec:
        src_path = fa_entry["path"]
        shutil.copy(
            os.path.join("out_pre_processed_fa", src_path + ".parquet.d", "chr21-nt.pdf"),
            os.path.join("fig", "chr21.d", fa_entry["name"] + ".chr21-nt.pdf")
        )
        shutil.copy(
            os.path.join("out_pre_processed_fa", src_path + ".parquet.d", "chrY-nt.pdf"),
            os.path.join("fig", "chrY.d", fa_entry["name"] + ".chrY-nt.pdf")
        )
