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
            PREPROCESSED_PATH = path("pre_processed_mtr", paste0(gtf_spec$path, ".mtr.parquet.d"))
        )
    if (gtf_spec$coordinate == "grch38"){
        all_gtf_spec_grch38 <- rbind(all_gtf_spec_grch38, this_gtf_spec)
    } else if (gtf_spec$coordinate == "chm13-t2t-v2"){
        all_gtf_spec_chm13_t2t_v2 <- rbind(all_gtf_spec_chm13_t2t_v2, this_gtf_spec)
    }
    rm(gtf_spec, this_gtf_spec)
}
mtr_jaccard <- arrow::read_parquet("summary_mtr_eq.parquet")

for (name in grep("_EQ", colnames(mtr_jaccard)[1:], value = TRUE)) {
    print(name)
    mtr_jaccard_grch38 <- mtr_jaccard %>%
        dplyr::inner_join(
            all_gtf_spec_grch38,
            by = c("F1" = "PREPROCESSED_PATH")
        ) %>%
        dplyr::rename(GTF_NAME1 = GTF_NAME) %>%
        dplyr::inner_join(
            all_gtf_spec_grch38,
            by = c("F2" = "PREPROCESSED_PATH")
        ) %>%
        dplyr::rename(GTF_NAME2 = GTF_NAME) %>%
        dplyr::mutate(
            DICE = 2 * .[[name]] / (LF1 + LF2)
        ) %>%
        dplyr::select(
            GTF_NAME1,
            GTF_NAME2,
            DICE
        )

    mtr_jaccard_grch38_wide <- mtr_jaccard_grch38 %>%
        tidyr::pivot_wider(
            names_from = GTF_NAME1,
            values_from = DICE
        )
    mtr_jaccard_grch38_wide_mtx <- mtr_jaccard_grch38_wide %>%
        as.data.frame()
    rownames(mtr_jaccard_grch38_wide_mtx) <- mtr_jaccard_grch38_wide_mtx$GTF_NAME2
    mtr_jaccard_grch38_wide_mtx$GTF_NAME2 <- NULL
    mtr_jaccard_grch38_wide_mtx <- as.matrix(mtr_jaccard_grch38_wide_mtx)

    pheatmap(
        mtr_jaccard_grch38_wide_mtx,
        color = c(
            grDevices::colorRampPalette(c('blue','white'))(100),
            grDevices::colorRampPalette(c('white','red'))(10)
        ),
        breaks = c(seq(0.0, 1.0, 0.01), seq(1.1, 10, 1.0)),
        display_numbers=TRUE,
        number_format="%.2f",
        filename = path("fig", sprintf("grch38_mtr_%s_jaccard.pdf", name))
    )
}
