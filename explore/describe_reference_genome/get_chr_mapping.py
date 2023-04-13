import json

import requests, sys
import pandas as pd

if __name__ == "__main__":
    server = "https://rest.ensembl.org"
    ext = "/info/assembly/homo_sapiens?"
    r = requests.get(
        server + ext,
        headers={"Content-Type": "application/json"},
        params={"synonyms": 1}
    )

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    with open("r.json", "w") as writer:
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
