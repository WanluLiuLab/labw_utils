import io
import json
import os
import sys

import pandas as pd
import requests

from labw_utils.commonutils.io.file_system import file_exists
from labw_utils.commonutils.io.safe_io import get_reader, get_writer


def get_ncbi_chromosome_spec():
    cache_path = os.path.join(CACHE_PATH, "ncbi_chromosome_cache.csv")
    if file_exists(cache_path):
        ncbi_df = pd.read_csv(cache_path)
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
            names="NCBIName	SequenceRole	AssignedMolecule	AssignedMolecule-Location	GenBankAccn	Relationship	RefSeqAccn	AssemblyUnit	SequenceLength	UCSC".split(
                "\t")
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
    cache_path = os.path.join(CACHE_PATH, "ensembl_chromosome_cache.json")
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

    retd = {}
    for contig in tlr:
        ucsc_name = None
        contig_names = [contig["name"]]
        for synonym in contig["synonyms"]:
            if synonym["dbname"] == "UCSC":
                ucsc_name = synonym["name"]
            else:
                contig_names.append(synonym["name"])
            if ucsc_name is not None:
                retd.update({contig_name:ucsc_name for contig_name in contig_names})
    return retd


def get_ucsc_chromosome_spec():
    cache_path = os.path.join(CACHE_PATH, "ucsc_chromosome_cache.txt")
    if file_exists(cache_path):
        with get_reader(cache_path) as reader:
            decoded = reader.read()
    else:
        r = requests.get(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chromAlias.txt"
        )
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.text
        with get_writer(cache_path) as writer:
            writer.write(decoded)
    retd = {}
    for line in decoded.splitlines():
        if line.startswith("#") or len(line) < 1:
            continue
        segments = line.split("\t")
        for i in range(1, len(segments)):
            if segments[i] == "":
                continue
            retd[segments[i]] = segments[0]

    return retd


if __name__ == "__main__":
    CACHE_PATH = "cache"
    os.makedirs(CACHE_PATH, exist_ok=True)
    ens_chr_spec = get_ensembl_chromosome_spec()
    ncbi_chr_spec = get_ncbi_chromosome_spec()
    ucsc_chr_spec = get_ucsc_chromosome_spec()
    final_spec = {}
    final_spec.update(ens_chr_spec)
    final_spec.update(ncbi_chr_spec)
    final_spec.update(ucsc_chr_spec)
    additional_spec = {} # To solve Ensembl bugs
    for k, v in final_spec.items():
        if k.startswith("HSCHR"):
            additional_spec["CHR_"+k] = v
    final_spec.update(additional_spec)
    with get_writer("final_chromosome_spec.json") as writer:
        json.dump(final_spec, writer, indent=4)
