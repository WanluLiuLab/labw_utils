"""
Preprocess -- Normalize chromosome name, keep transcript and exon feature only.
"""
import gc
import glob
import json
import os.path

from labw_utils.bioutils.datastructure.feature_view import FeatureView
from labw_utils.commonutils.lwio.safe_io import get_reader


if __name__ == "__main__":
    with get_reader("final_chromosome_spec.json") as reader:
        final_chromosome_spec = json.load(reader)

    for fn in glob.glob(os.path.join("gtf", "*.gtf")):
        _ = FeatureView.from_gtf(fn)
        gc.collect()
