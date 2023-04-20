library(rjson)
library(tidyverse)
library(fs)

gtf_spec_json <- rjson::fromJSON(file = "gtf_spec.json")
all_gtf_spec <- data.frame()
for (gtf_spec in gtf_spec_json) {
    all_gtf_spec <- rbind(
        all_gtf_spec,
        data.frame(
            GTF_NAME = gtf_spec$name,
            PREPROCESSED_PATH = path("pre_processed_gtf", paste0(gtf_spec$path))
        )
    )
    rm(gtf_spec)
}

gtf_chr_isoform_cnt <- readr::read_tsv(
    "gtf_chr_cnt.tsv",
    col_names = c("PREPROCESSED_PATH", "CHR", "NUM")
)
gtf_chr_isoform_cnt_merged <- gtf_chr_isoform_cnt %>%
    dplyr::inner_join(
        all_gtf_spec,
        by = "PREPROCESSED_PATH"
    )
gtf_chr_isoform_cnt_merged_chrUn <- gtf_chr_isoform_cnt_merged %>%
    dplyr::transmute(
        CHR_REPLACED = str_replace_all(
            .$CHR,
            ".*_.*",
            "chrUn"
        ),
        NUM,
        GTF_NAME,
        CHR_ISNORMAL = str_count(
            .$CHR,
            "_"
        ) == 0
    ) %>%
    dplyr::group_by(
        GTF_NAME, CHR_REPLACED
    ) %>%
    dplyr::mutate(
        NUM = sum(NUM)
    ) %>%
    unique()

gtf_chr_isoform_cnt_merged_chrUn_sum <- gtf_chr_isoform_cnt_merged_chrUn %>%
    dplyr::group_by(
        GTF_NAME, CHR_ISNORMAL
    ) %>%
    dplyr::mutate(
        NUM = sum(NUM)
    ) %>%
    dplyr::select(NUM, GTF_NAME, CHR_ISNORMAL) %>%
    unique()

g <- ggplot(gtf_chr_isoform_cnt_merged_chrUn_sum) +
    geom_bar(
        aes(
            x = NUM,
            y = GTF_NAME,
            fill = CHR_ISNORMAL
        ),
        stat = "identity"
    ) +
    scale_fill_manual(
        "Is chr[1-23XYM]?",
        limits = c(FALSE, TRUE),
        values = c("red", "black")
    ) +
    theme_bw() +
    ggtitle("N. Isoforms inside each GTF")
ggsave(path("fig", "gtf-chr-len.svg"), g, width = 10, height = 5)

gtf_chr_usage <- gtf_chr_isoform_cnt_merged %>%
    dplyr::select(GTF_NAME, CHR) %>%
    dplyr::mutate(CHR_PRESENCE = TRUE) %>%
    tidyr::pivot_wider(
        names_from=GTF_NAME,
        values_from=CHR_PRESENCE,
        id_cols=CHR,
        values_fill = FALSE
    )
gtf_chr_usage_df <- as.data.frame(gtf_chr_usage)
rownames(gtf_chr_usage_df) <- gtf_chr_usage_df$CHR
gtf_chr_usage_df$CHR <- NULL

g <- ComplexUpset::upset(gtf_chr_usage_df, colnames(gtf_chr_usage_df))
ggsave(path("fig", "gtf-chr-upset.svg"), g, width=15, height=8)
