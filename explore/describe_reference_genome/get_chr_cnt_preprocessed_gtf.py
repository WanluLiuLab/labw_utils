import glob
import multiprocessing
import os
from collections import defaultdict

from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.feature import FeatureType


def count_chr(
        src_gtf_path: str
):
    retd = defaultdict(lambda: 0)
    for record in GtfIterator(src_gtf_path):
        if record.parsed_feature == FeatureType.TRANSCRIPT:
            retd[record.seqname] += 1
    for seqname, cnt in retd.items():
        print("\t".join((src_gtf_path, seqname, str(cnt))))


if __name__ == "__main__":
    ppool = []
    for fn in glob.glob(os.path.join("pre_processed_gtf", "*.gtf")):
        ppool.append(multiprocessing.Process(target=count_chr, args=(fn,)))
        ppool[-1].start()
    for p in ppool:
        p.join()
