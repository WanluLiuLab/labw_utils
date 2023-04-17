"""
Get sample attributres from MTRs.

TODO
"""
import glob
import json
import os.path

import pandas as pd

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer

if __name__ == "__main__":
    retd = {}
    for fn in tqdm(list(glob.glob(os.path.join("pre_processed_mtr", "*.parquet")))):
        df = pd.read_parquet(fn)
        retd[fn] = list(
            json.loads(attr) for attr in df.sample(min(1000, len(df)))["ATTRS"]
        )
    with get_writer("final_mtr_attr_sample.json") as writer:
        json.dump(retd, writer, indent=4)
