library(ComplexUpset)
library(tidyverse)
library(rjson)
library(fs)
fa_spec_json <- rjson::fromJSON(file = "fa_spec.json")
all_fa_spec <- data.frame()
for (fa_spec in fa_spec_json) {
    all_fa_spec <- rbind(
        all_fa_spec,
        data.frame(
            FA_NAME = fa_spec$name,
            METADATA_PATH = path("out_pre_processed_fa", paste0(fa_spec$path, ".json")),
            FAI_PATH = path("pre_processed_fa", paste0(fa_spec$path, ".fai"))
        )
    )
    rm(fa_spec)
}

all_chr_spec <- data.frame()


for (fn in Sys.glob(path("out_pre_processed_fa", "*.json"))) {
    message(fn)
    fc <- rjson::fromJSON(file = fn)
    for (chr_spec in fc$FASTA_CHRS) {
        if (is.null(chr_spec$TYPE)) {
            chr_type <- "NULL"
        }
        else if (str_count(chr_spec$TYPE$toplevel, "Analysis Set") > 0) {
            chr_type <- paste(chr_spec$TYPE$toplevel, chr_spec$TYPE$details$TYPE, sep = ":")
        } else {
            chr_type <- chr_spec$TYPE$toplevel
        }
        all_chr_spec <- rbind(
            all_chr_spec,
            data.frame(
                NAME = chr_spec$NAME,
                LEN = chr_spec$LEN,
                TYPE = chr_type,
                FN = fn
            )
        )
        rm(chr_spec, chr_type)
    }
    rm(fn, fc)
}

all_chr_spec_joint <- all_chr_spec %>%
    dplyr::inner_join(
        all_fa_spec,
        by = c("FN" = "METADATA_PATH")
    )
all_chr_spec_uniq <- all_chr_spec %>%
    dplyr::select(NAME, TYPE) %>%
    unique()

ggplot(all_chr_spec_joint) +
    geom_bar(
        aes(
            y = FA_NAME,
            fill = TYPE
        )
    ) +
    theme_bw()
ggsave(path("fig", "fasta-chr-type.svg"), width = 10, height = 5)


ggplot(all_chr_spec_joint) +
    geom_bar(
        aes(
            x = LEN,
            y = FA_NAME,
            fill = TYPE
        ),
        stat = "identity"
    ) +
    scale_x_continuous(
        labels = scales::label_number(
            scale_cut = scales::cut_long_scale()
        )
    ) +
    theme_bw()
ggsave(path("fig", "fasta-chr-type-len.svg"), width = 10, height = 5)


all_fai <- data.frame()

for (fn in Sys.glob(path("pre_processed_fa", "*.fai"))) {
    message(fn)
    fai <- readr::read_tsv(fn, col_names = c("CHR", "LEN", "OFFSET", "LINE_BLEN", "LINE_LEN")) %>%
        dplyr::transmute(
            CHR,
            FAI_PATH = fn
        )
    all_fai <- rbind(all_fai, fai)
    rm(fai, fn)
}

all_fai_merged <- all_fai %>%
    dplyr::inner_join(
        all_fa_spec,
        by = "FAI_PATH"
    ) %>%
    dplyr::select(!FAI_PATH) %>%
    dplyr::mutate(CHR_PRESENCE = TRUE) %>%
    tidyr::pivot_wider(
        names_from = FA_NAME,
        values_from = CHR_PRESENCE,
        id_cols = CHR,
        values_fill = FALSE
    ) %>%
    dplyr::inner_join(
        all_chr_spec_uniq,
        by = c("CHR" = "NAME")
    )
all_fai_merged_df <- as.data.frame(all_fai_merged)
rownames(all_fai_merged_df) <- all_fai_merged_df$CHR
all_fai_merged_df$CHR <- NULL

g <- ComplexUpset::upset(
    data=all_fai_merged_df,
    intersect= setdiff(colnames(all_fai_merged_df), "TYPE"),
    base_annotations = list(
        'Intersection size' = intersection_size(
            mapping = aes(fill = TYPE)
        )),
    set_sizes = (
        upset_set_size(
            geom = geom_bar(
                aes(fill = TYPE, group = TYPE),
                width = 0.8
            )
        )
    ),
    guides = 'over'
)

ggsave(path("fig", "fasta-chr-upset.svg"), g, width = 12, height = 10)
