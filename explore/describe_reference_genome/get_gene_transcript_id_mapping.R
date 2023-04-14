library("tidyverse")
library("arrow")
library("biomaRt")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

getBM(
    attributes = c("ensembl_gene_id_version", "ensembl_transcript_id_version"),
    mart = ensembl
) %>%
    arrow::write_parquet("ens_gene_transcript_id_map.parquet")

getBM(
    attributes = c("ensembl_gene_id_version", "hgnc_id", "hgnc_symbol"),
    mart = ensembl
) %>%
    arrow::write_parquet("ens_hgnc_gene_map.parquet")

getBM(
    attributes = c("ensembl_transcript_id_version", "hgnc_trans_name"),
    mart = ensembl
) %>%
    arrow::write_parquet("ens_hgnc_transcript_map.parquet")

getBM(
    attributes = c(
        "ensembl_transcript_id_version",
        "refseq_mrna",
        "refseq_ncrna",
        "refseq_peptide"
    ),
    mart = ensembl
) %>%
    arrow::write_parquet("ens_refseq_transcript_map.parquet")
# ENST00000641006.1
getBM(
    attributes = c(
        "ensembl_transcript_id_version",
        "refseq_mrna_predicted",
        "refseq_ncrna_predicted",
        "refseq_peptide_predicted"
    ),
    mart = ensembl
) %>%
    arrow::write_parquet("ens_refseq_predicted_transcript_map.parquet")


