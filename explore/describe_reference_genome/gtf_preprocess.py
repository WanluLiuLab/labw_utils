"""
Preprocess -- Normalize chromosome name, keep transcript and exon feature only.
"""
import glob
import json
import multiprocessing
import os.path

from labw_utils.bioutils.accession_matcher import infer_accession_type
from labw_utils.bioutils.parser.gtf import GtfIterator, GtfIteratorWriter
from labw_utils.bioutils.record.feature import FeatureType
from labw_utils.commonutils.io.safe_io import get_reader, get_writer


def preprocess(
        src_path: str,
        dst_path: str
):
    not_recognized_chr = set()
    n_pass = 0
    n_filtered_no_recognized_feature = 0
    n_filtered_no_recognized_chr = 0
    with GtfIterator(src_path, show_tqdm=True) as reader:
        with GtfIteratorWriter(dst_path) as writer:
            for record in reader:
                if record.parsed_feature not in (
                        FeatureType.EXON,
                        FeatureType.TRANSCRIPT
                ):
                    n_filtered_no_recognized_feature += 1
                    continue
                chr_name = record.seqname
                chr_name_type = infer_accession_type(chr_name)
                if chr_name_type is not None and chr_name_type.toplevel.startswith("Analysis Set"):
                    n_pass += 1
                    writer.write(record)
                    continue
                converted_name = final_chromosome_spec.get(chr_name, None)
                if converted_name is not None:
                    converted_name_type = infer_accession_type(converted_name)
                    if converted_name_type is not None and converted_name_type.toplevel.startswith("Analysis Set"):
                        record = record.update(seqname=converted_name)
                        n_pass += 1
                        writer.write(record)
                        continue
                not_recognized_chr.add(chr_name)
                n_filtered_no_recognized_chr += 1
    with get_writer(dst_path + ".norecognized.txt") as writer:
        writer.write("\n".join(list(not_recognized_chr)))
    print(f"{src_path} -> {dst_path}: SUCCESS {n_pass} NO_RECOGNIZED_CHR {n_filtered_no_recognized_chr}")


if __name__ == "__main__":
    with get_reader("final_chromosome_spec.json") as reader:
        final_chromosome_spec = json.load(reader)
    os.makedirs("pre_processed_gtf", exist_ok=True)

    ppool = []
    for fn in glob.glob(os.path.join("gtf", "*.gtf")):
        dst_fn = os.path.join("pre_processed_gtf", os.path.basename(fn))
        ppool.append(multiprocessing.Process(target=preprocess, args=(fn, dst_fn)))
        ppool[-1].start()
    for p in ppool:
        p.join()
