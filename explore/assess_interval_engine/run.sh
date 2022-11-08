#!/usr/bin/env bash
# shellcheck disable=SC1091

# if [ ! -e venv ]; then
#     python3 -m venv venv
#     source venv/bin/activate
#     pip install -r requirements.txt
# fi

python3 ./extract_from_bam.py /home/yuzj/Documents/cpptetgs_experimental/test_data/real/cDNA_sorted.bam bam_regions.tsv
bash ./extract_gtf.sh /home/yuzj/Documents/cpptetgs_experimental/test_data/real/Homo_sapiens.GRCh38.105.gtf gtf_regions.tsv

pv bam_regions.tsv | grep -e "^chr" -e "^1[[:space:]]" > bam_regions_chr1.tsv
pv gtf_regions.tsv | grep -e "^chr" -e "^1[[:space:]]" > gtf_regions_chr1.tsv
