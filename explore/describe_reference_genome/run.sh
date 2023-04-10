#!/usr/bin/env bash

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

Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/ncbi_analysis_set.parquet \
    --output out/ncbi_analysis_set-plots/

mkdir -p out/gencode_toplevel_comprehensive-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f GRCh38.p13.genome.fa \
    -g gencode.v43.chr_patch_hapl_scaff.annotation.gtf \
    -o out/gencode_toplevel_comprehensive
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/gencode_toplevel_comprehensive.parquet \
    --output out/gencode_toplevel_comprehensive-plots/

mkdir -p out/ncbi_refseq-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f GCF_000001405.40_GRCh38.p14_genomic.fna \
    -g GCF_000001405.40_GRCh38.p14_genomic.gtf \
    -o out/ncbi_refseq
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/ncbi_refseq.parquet \
    --output out/ncbi_refseq-plots/

mkdir -p out/ncbi_genebank-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f GCA_000001405.29_GRCh38.p14_genomic.fna \
    -g GCA_000001405.29_GRCh38.p14_genomic.gtf \
    -o out/ncbi_genebank
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/ncbi_genebank.parquet \
    --output out/ncbi_genebank-plots/

mkdir -p out/gencode_primary_assembly_comprehensive-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f GRCh38.primary_assembly.genome.fa \
    -g gencode.v43.primary_assembly.annotation.gtf \
    -o out/gencode_primary_assembly_comprehensive
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/gencode_primary_assembly_comprehensive.parquet \
    --output out/gencode_primary_assembly_comprehensive-plots/

mkdir -p out/gencode_toplevel_basic-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f GRCh38.p13.genome.fa \
    -g gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtf \
    -o out/gencode_toplevel_basic
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/gencode_toplevel_basic.parquet \
    --output out/gencode_toplevel_basic-plots/

mkdir -p out/gencode_primary_assembly_basic-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f GRCh38.primary_assembly.genome.fa \
    -g gencode.v43.primary_assembly.basic.annotation.gtf \
    -o out/gencode_primary_assembly_basic
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/gencode_primary_assembly_basic.parquet \
    --output out/gencode_primary_assembly_basic-plots/

mkdir -p out/ucsc_analysis_set-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f hg38.analysisSet.fa \
    -g hg38.ncbiRefSeq.gtf \
    -o out/ucsc_analysis_set
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/ucsc_analysis_set.parquet \
    --output out/ucsc_analysis_set-plots/

mkdir -p out/usuc_defaults-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f hg38.fa \
    -g hg38.ncbiRefSeq.gtf \
    -o out/usuc_defaults
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/ucsc_analysis_set.parquet \
    --output out/usuc_defaults-plots/

mkdir -p out/usuc_masked-plots
python -m labw_utils.bioutils describe_fasta_gtf_by_binning \
    -f hg38.masked.fa \
    -g hg38.ncbiRefSeq.gtf \
    -o out/usuc_masked
Rscript plot_describe_fasta_gtf_by_binning.R \
    --input out/ucsc_analysis_set.parquet \
    --output out/usuc_masked-plots/
