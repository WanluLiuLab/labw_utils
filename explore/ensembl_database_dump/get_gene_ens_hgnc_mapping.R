library("tidyverse")
library("arrow")
library("biomaRt")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

attrs <- listAttributes(ensembl) %>% as_tibble()

getBM(
    attributes = c("ensembl_gene_id", "version", "hgnc_id", "hgnc_symbol"),
    mart = ensembl
) %>%
    arrow::write_parquet("ens_hgnc_gene_map.parquet")
