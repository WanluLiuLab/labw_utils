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

mkdir -p pre_processed_gtf_bedtools_sorted
for fn in pre_processed_gtf/*.gtf; do
    bedtools sort -i "${fn}" > pre_processed_gtf_bedtools_sorted/"$(basename "${fn}")" &
done
wait

printf "F1\tF2\tJaccard\n" > bedtools_jaccard.tsv
for afn in pre_processed_gtf_bedtools_sorted/*.gtf; do
    for bfn in pre_processed_gtf_bedtools_sorted/*.gtf; do
        printf "%s\t%s\t" "${afn}" "${bfn}" > bedtools_jaccard.tsv
        bedtools jaccard -s -a "${afn}" -b "${bfn}" | tail -n 1 | cut -f 3 > bedtools_jaccard.tsv
    done
done
