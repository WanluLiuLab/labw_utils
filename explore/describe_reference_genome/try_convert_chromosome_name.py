import json

from labw_utils.bioutils.accession_matcher import infer_accession_type
from labw_utils.commonutils.io.file_system import file_exists
from labw_utils.commonutils.io.safe_io import get_reader
from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

setup_basic_logger()

_lh = get_logger()

if __name__ == "__main__":
    with get_reader("fa_spec.json") as reader:
        fa_spec = json.load(reader)
    with get_reader("final_chromosome_spec.json") as reader:
        final_chromosome_spec = json.load(reader)
    for fa_entry in fa_spec:
        metadata_json_path = "out_" + fa_entry["path"] + ".json"
        if file_exists(metadata_json_path):
            with get_reader(metadata_json_path) as reader:
                metadata = json.load(reader)
                is_analysis_set = 0
                unknown = 0
                for chr_entry in metadata["FASTA_CHRS"]:
                    chr_name = chr_entry["NAME"]
                    chr_name_type = infer_accession_type(chr_name)
                    if chr_name_type is not None and chr_name_type.toplevel.startswith("Analysis Set"):
                        is_analysis_set += 1
                        continue
                    converted_name = final_chromosome_spec.get(chr_name, None)
                    if converted_name is not None:
                        converted_name_type = infer_accession_type(converted_name)
                        if converted_name_type is not None and converted_name_type.toplevel.startswith("Analysis Set"):
                            is_analysis_set += 1
                            continue

                    unknown += 1
                    print(fa_entry["name"], chr_name)
        else:
            _lh.warning("Name %s skipped", fa_entry["name"])
