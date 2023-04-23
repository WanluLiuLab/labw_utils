#!/usr/bin/env bash
# TODO: Illumina iGenome https://support.illumina.com/sequencing/sequencing_software/igenome.html

export PYTHONPATH=/home/yuzj/Documents/labw_utils/src

axel https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc_complete_set.json -o id_map/hgnc_complete_set.json
axel https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/withdrawn.json -o id_map/hgnc_withdrawn.json
python get_chr_mapping.py
Rscript get_gene_transcript_id_mapping.R

mkdir -p out out_fa out_gtf out_pre_processed_fa fig out/ncbi_analysis_set-plots fig/chr21.d fig/chrY.d

for fn in fa/*.fa fa/*.fna; do
    python -m labw_utils.bioutils describe_fasta_by_binning \
        -f "${fn}" \
        -o out_"${fn}" \
        --metadata_only &
done
wait
python fa_preprocess.py
for fn in pre_processed_fa/Homo_sapiens.GRCh38.dna*.fa ; do # pre_processed_fa/*.fna
    {
        samtools faidx "${fn}"
        python -m labw_utils.bioutils describe_fasta_by_binning \
            -f "${fn}" \
            -o out_"${fn}"
    } &
done
wait

for fn in out_pre_processed_fa/*.parquet; do
    Rscript plot_fasta_base_usage.R \
        --input "${fn}" \
        --output "${fn}.d"
done
python fa_copy_chr21.py
Rscript plot_fasta_chromosome_type.R

python gtf_preprocess.py
python gtf_to_mtr.py
python get_sample_mtr_attr.py
python get_chr_cnt_preprocessed_gtf.py >gtf_chr_cnt.tsv
