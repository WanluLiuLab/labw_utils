library(tidyverse)
library(rjson)
library(fs)
library(corrplot)
library(pheatmap)

gtf_spec_json <- rjson::fromJSON(file = "gtf_spec.json")
all_gtf_spec_grch38 <- data.frame()
all_gtf_spec_chm13_t2t_v2 <- data.frame()

for (gtf_spec in gtf_spec_json) {
    this_gtf_spec <- data.frame(
            GTF_NAME = gtf_spec$name,
            PREPROCESSED_PATH = path("pre_processed_gtf_bedtools_sorted", paste0(gtf_spec$path))
        )
    if (gtf_spec$coordinate == "grch38"){
        all_gtf_spec_grch38 <- rbind(all_gtf_spec_grch38, this_gtf_spec)
    } else if (gtf_spec$coordinate == "chm13-t2t-v2"){
        all_gtf_spec_chm13_t2t_v2 <- rbind(all_gtf_spec_chm13_t2t_v2, this_gtf_spec)
    }
    rm(gtf_spec, this_gtf_spec)
}
bedtools_jaccard <- readr::read_tsv("bedtools_jaccard.tsv")

bedtools_jaccard_grch38 <- bedtools_jaccard %>%
    dplyr::inner_join(
        all_gtf_spec_grch38,
        by=c("F1"="PREPROCESSED_PATH")
    ) %>%
    dplyr::rename(GTF_NAME1=GTF_NAME) %>%
    dplyr::inner_join(
        all_gtf_spec_grch38,
        by=c("F2"="PREPROCESSED_PATH")
    ) %>%
    dplyr::rename(GTF_NAME2=GTF_NAME) %>%
    dplyr::select(
        GTF_NAME1,
        GTF_NAME2,
        Jaccard
    )
bedtools_jaccard_chm13_t2t_v2 <- bedtools_jaccard %>%
    dplyr::inner_join(
        all_gtf_spec_chm13_t2t_v2,
        by=c("F1"="PREPROCESSED_PATH")
    ) %>%
    dplyr::rename(GTF_NAME1=GTF_NAME) %>%
    dplyr::inner_join(
        all_gtf_spec_chm13_t2t_v2,
        by=c("F2"="PREPROCESSED_PATH")
    ) %>%
    dplyr::rename(GTF_NAME2=GTF_NAME) %>%
    dplyr::select(
        GTF_NAME1,
        GTF_NAME2,
        Jaccard
    )
bedtools_jaccard_grch38_wide <- bedtools_jaccard_grch38 %>%
    tidyr::pivot_wider(
        names_from = GTF_NAME1,
        values_from = Jaccard
    )
bedtools_jaccard_grch38_wide_mtx <- bedtools_jaccard_grch38_wide %>%
    as.data.frame()
rownames(bedtools_jaccard_grch38_wide_mtx) <- bedtools_jaccard_grch38_wide_mtx$GTF_NAME2
bedtools_jaccard_grch38_wide_mtx$GTF_NAME2 <- NULL
bedtools_jaccard_grch38_wide_mtx <- as.matrix(bedtools_jaccard_grch38_wide_mtx)
pheatmap(
    bedtools_jaccard_grch38_wide_mtx,
    filename = path("fig", "grch38_gtf_jaccard.pdf"),
    breaks=seq(0, 1, 0.01)
)

bedtools_jaccard_chm13_t2t_v2_wide <- bedtools_jaccard_chm13_t2t_v2 %>%
    tidyr::pivot_wider(
        names_from = GTF_NAME1,
        values_from = Jaccard
    )
bedtools_jaccard_chm13_t2t_v2_wide_mtx <- bedtools_jaccard_chm13_t2t_v2_wide %>%
    as.data.frame()
rownames(bedtools_jaccard_chm13_t2t_v2_wide_mtx) <- bedtools_jaccard_chm13_t2t_v2_wide_mtx$GTF_NAME2
bedtools_jaccard_chm13_t2t_v2_wide_mtx$GTF_NAME2 <- NULL
bedtools_jaccard_chm13_t2t_v2_wide_mtx <- as.matrix(bedtools_jaccard_chm13_t2t_v2_wide_mtx)
pheatmap(
    bedtools_jaccard_chm13_t2t_v2_wide_mtx,
    filename = path("fig", "chm13_t2t_v2_gtf_jaccard.pdf"),
    breaks=seq(0, 1, 0.01)
)
