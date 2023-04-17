#!/usr/bin/env bash
# TODO: Illumina iGenome https://support.illumina.com/sequencing/sequencing_software/igenome.html

export PYTHONPATH=/home/yuzj/Documents/labw_utils/src

mkdir -p out out_fa out_gtf
# make

mkdir -p out/ncbi_analysis_set-plots
for fn in fa/*.fa fa/*.fna; do
    python -m labw_utils.bioutils describe_fasta_by_binning \
        -f "${fn}" \
        -o out_"${fn}" \
        --metadata_only
done

cat pre_processed_gtf/*.gtf | cut -f 2 | sort | uniq >all_sources.txt

axel https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc_complete_set.json
axel https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/withdrawn.json -o hgnc_withdrawn.json
