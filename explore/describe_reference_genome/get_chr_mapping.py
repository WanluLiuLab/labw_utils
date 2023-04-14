import io
import json

import requests, sys
import pandas as pd

from labw_utils.commonutils.io.file_system import file_exists
from labw_utils.commonutils.io.safe_io import get_reader, get_writer


def get_ncbi_chromosome_spec():
    cache_path = "ncbi_chromosome_cache.csv"
    if file_exists(cache_path):
        ncbi_df =  pd.read_csv(cache_path)
    else:
        url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/' \
              'GCF_000001405.40_GRCh38.p14_assembly_report.txt'
        r = requests.get(url)
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        ncbi_df = pd.read_table(
            io.StringIO(r.text),
            comment="#",
            names="NCBIName	SequenceRole	AssignedMolecule	AssignedMolecule-Location	GenBankAccn	Relationship	RefSeqAccn	AssemblyUnit	SequenceLength	UCSC".split("\t")
        )
        ncbi_df.to_csv(cache_path, index=False)
    retd = {}

    synonym_dbnames_no_ucsc = ["NCBIName", "GenBankAccn", "RefSeqAccn"]
    for df_tuple in ncbi_df.itertuples():
        ucsc_value = df_tuple.UCSC
        if ucsc_value is not None:
            for synonym_dbname_no_ucsc in synonym_dbnames_no_ucsc:
                db_value = df_tuple.__getattribute__(synonym_dbname_no_ucsc)
                if db_value is not None:
                    retd[db_value] = ucsc_value
    return retd




def get_ensembl_chromosome_spec():
    cache_path = "ensembl_chromosome_cache.json"
    if file_exists(cache_path):
        with get_reader(cache_path) as reader:
            decoded = json.load(reader)
    else:
        r = requests.get(
            "https://rest.ensembl.org/info/assembly/homo_sapiens?",
            headers={"Content-Type": "application/json"},
            params={"synonyms": 1}
        )
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        with get_writer(cache_path) as writer:
            json.dump(decoded, writer, indent=4)

    tlr = decoded["top_level_region"]
    synonym_dbnames = set()
    for contig in tlr:
        for synonym in contig["synonyms"]:
            synonym_dbnames.add(synonym["dbname"])
    synonym_dbnames = list(synonym_dbnames)
    final_table = {
        "ens_name": [],
        **{synonym_dbname: [] for synonym_dbname in synonym_dbnames}
    }
    for contig in tlr:
        this_synonym_dbnames = list(synonym_dbnames)
        for synonym in contig["synonyms"]:
            if synonym["name"] != "na":
                final_table[synonym["dbname"]].append(synonym["name"])
                this_synonym_dbnames.remove(synonym["dbname"])
        for na_dbname in this_synonym_dbnames:
            final_table[na_dbname].append(None)

        final_table["ens_name"].append(contig["name"])
    final_table_df = pd.DataFrame(final_table)
    retd = {}

    synonym_dbnames_no_ucsc = list(synonym_dbnames)
    synonym_dbnames_no_ucsc.remove("UCSC")
    synonym_dbnames_no_ucsc.append("ens_name")
    for df_tuple in final_table_df.itertuples():
        ucsc_value = df_tuple.UCSC
        if ucsc_value is not None:
            for synonym_dbname_no_ucsc in synonym_dbnames_no_ucsc:
                db_value = df_tuple.__getattribute__(synonym_dbname_no_ucsc)
                if db_value is not None:
                    retd[db_value] = ucsc_value
    return retd


if __name__ == "__main__":
    ens_chr_spec = get_ensembl_chromosome_spec()
    ncbi_chr_spec = get_ncbi_chromosome_spec()
    final_spec = {}
    final_spec.update(ens_chr_spec)
    final_spec.update(ncbi_chr_spec)
    with get_writer("final_chromosome_spec.json") as writer:
        json.dump(final_spec, writer, indent=4)
