import json
import multiprocessing
import os

from labw_utils.bioutils.accession_matcher import infer_accession_type
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.lwio.safe_io import get_reader
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Mapping

setup_basic_logger()

_lh = get_logger()


def convert(
        final_chromosome_spec: Mapping[str, str],
        fa_entry: Mapping[str, str]
):
    src_path = os.path.join("fa", fa_entry["path"])
    dst_path = os.path.join("pre_processed_fa", fa_entry["path"])
    fav = FastaViewFactory(src_path, read_into_memory=False)
    with FastaWriter(dst_path, split_at=80) as writer:
        for chr_name in tqdm(fav.chr_names):
            chr_name_type = infer_accession_type(chr_name)
            if chr_name_type is not None and chr_name_type.toplevel.startswith("Analysis Set"):
                converted_name = chr_name
            else:
                converted_name = final_chromosome_spec.get(chr_name, None)
                if converted_name is not None:
                    converted_name_type = infer_accession_type(converted_name)
                    if converted_name_type is not None and converted_name_type.toplevel.startswith(
                            "Analysis Set"
                    ):
                        pass
                    else:
                        continue
                else:
                    continue
            writer.write(
                FastaRecord(
                    seq_id=converted_name,
                    sequence=fav.sequence(chr_name)
                )
            )


if __name__ == "__main__":
    os.makedirs("pre_processed_fa", exist_ok=True)
    ppool = []
    with get_reader("fa_spec.json") as reader:
        fa_spec = json.load(reader)
    with get_reader("final_chromosome_spec.json") as reader:
        _final_chromosome_spec = json.load(reader)
    for _fa_entry in fa_spec:
        ppool.append(multiprocessing.Process(
            target=convert,
            args=(_final_chromosome_spec, _fa_entry)
        ))
        ppool[-1].start()

    for p in ppool:
        p.join()
