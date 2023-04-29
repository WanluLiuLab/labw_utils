library("tidyverse")
library("arrow")
library("biomaRt")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

attrs <- listAttributes(ensembl) %>% as_tibble()

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

for (refseq_name in c(
    "refseq_mrna",
    "refseq_ncrna",
    "refseq_mrna_predicted",
    "refseq_ncrna_predicted",
    "refseq_peptide_predicted",
    "refseq_peptide"
)) {
    getBM(
        attributes = c(
            "ensembl_transcript_id_version",
            refseq_name
        ),
        mart = ensembl
    ) %>%
        arrow::write_parquet(sprintf("ens_%s_transcript_map.parquet", refseq_name))
}
