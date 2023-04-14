"""
Get sample attributres from MTRs.

TODO
"""
import glob
import os.path

import pandas as pd

if __name__ == "__main__":

    dfs = []
    for fn in glob.glob(os.path.join("pre_processed_mtr", "*.parquet")):
        df = pd.read_parquet(fn).sample(1000).loc[:, "ATTRS", "FPATH"]
        dfs.append(df)
    pd.concat(dfs).to_("final_mtr.parquet", index=False)
