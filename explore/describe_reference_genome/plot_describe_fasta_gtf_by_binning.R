library(argparser)

p <- arg_parser("")
p <- add_argument(p, "--input", help = "input Parquet file")
p <- add_argument(p, "--output", help = "Output Directory")
argv <- parse_args(p)

library(tidyverse)
library(arrow)
library(parallel)

cl <- parallel::makeCluster(40)
df <- arrow::read_parquet(argv$input) %>%
    tibble::as_tibble() %>%
    replace(is.na(.), 0)
df_cat <- df %>%
    dplyr::select(!c(
        gtf_intervals_pos,
        gtf_intervals_neg,
        gtf_intervals_strandless,
        sdi,
        gc
    ))
df_gtf_intervals <- df %>%
    dplyr::select(
        start,
        chr_name,
        gtf_intervals_pos,
        gtf_intervals_neg,
        gtf_intervals_strandless
    )
df_diversity <- df %>%
    dplyr::select(
        start,
        chr_name,
        sdi,
        gc
    )


normal_nts <- c("A", "C", "G", "T")
sm_nts <- tolower(normal_nts)
hm_nts <- "N"
sm_hm_nts <- tolower(hm_nts)
ambig_nts <- c(
    "R", #	AG	puRine
    "Y", #	CT	pYrimidine
    "K", #	GT	Keto
    "M", #	AC	aMino
    "S", #	GC	Strong
    "W", #	AT	Weak
    "B", #	CGT	Not A
    "D", #	AGT	Not C
    "H", #	ACT	Not G
    "V" #	ACG	Not T
)
sm_ambig_nts <- tolower(ambig_nts)
other_nts <- setdiff(
    colnames(df_cat),
    c(
        normal_nts,
        sm_nts,
        hm_nts,
        sm_hm_nts,
        ambig_nts,
        sm_ambig_nts,
        "chr_name",
        "start"
    )
)

nt_scale_limits <- c(
    normal_nts,
    sm_nts,
    hm_nts,
    sm_hm_nts,
    ambig_nts,
    sm_ambig_nts
)
nt_scale_colors <- c(
    c("green", "red", "yellow", "cyan"),
    c("green4", "red4", "yellow4", "cyan4"),
    "black",
    grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(5, "Greys")
    )(length(c(sm_hm_nts, ambig_nts, sm_ambig_nts)))
)

nt_cat_scale_limits <- c(
    "NT_NORMAL",
    "NT_SOFT_MASKED",
    "NT_HARD_MASKED",
    "NT_SOFT_HARD_MASKED",
    "NT_AMBIG",
    "NT_AMBIG_SOFT_MASKED",
    "OTHERS"
)
nt_cat_scale_colors <- RColorBrewer::brewer.pal(7, "Set1")

parallel::clusterExport(cl, varlist = ls())

plot_chr <- function(this_chr_name) {
    library(tidyverse)
    this_df_cat <- df_cat %>%
        dplyr::filter(chr_name == this_chr_name) %>%
        dplyr::select(!chr_name)
    this_df_cat_long <- this_df_cat %>%
        tidyr::gather(
            key = "Nt",
            value = "value",
            -start
        )
    this_df_cat_long_1 <- this_df_cat %>%
        dplyr::transmute(
            NT_NORMAL = rowSums(dplyr::select(this_df_cat, dplyr::any_of(normal_nts))),
            NT_SOFT_MASKED = rowSums(dplyr::select(this_df_cat, dplyr::any_of(sm_nts))),
            NT_HARD_MASKED = rowSums(dplyr::select(this_df_cat, dplyr::any_of(hm_nts))),
            NT_SOFT_HARD_MASKED = rowSums(dplyr::select(this_df_cat, dplyr::any_of(sm_nts))),
            NT_AMBIG = rowSums(dplyr::select(this_df_cat, dplyr::any_of(ambig_nts))),
            NT_AMBIG_SOFT_MASKED = rowSums(dplyr::select(this_df_cat, dplyr::any_of(sm_ambig_nts))),
            OTHERS = rowSums(dplyr::select(this_df_cat, dplyr::any_of(other_nts))),
            start = start
        ) %>%
        tidyr::gather(
            key = "Nt",
            value = "value",
            -start
        )
    this_df_gtf_intervals <- df_gtf_intervals %>%
        dplyr::filter(chr_name == this_chr_name) %>%
        dplyr::select(!chr_name) %>%
        tidyr::gather(
            key = "IntervalType",
            value = "value",
            -start
        )
    this_df_diversity <- df_diversity %>%
        dplyr::filter(chr_name == this_chr_name) %>%
        dplyr::select(!chr_name) %>%
        tidyr::gather(
            key = "DiversityIndexName",
            value = "value",
            -start
        )

    message(sprintf("Plotting %s: Plotting", this_chr_name))
    g1 <- ggplot(this_df_cat_long) +
        geom_area(
            aes(x = start, y = value, fill = Nt)
        ) +
        scale_x_continuous("POS", labels = scales::label_number()) +
        scale_fill_manual(
            limits = nt_scale_limits,
            values = nt_scale_colors
        ) +
        theme_bw() +
        ggtitle(sprintf("Base distribution in %s", this_chr_name))

    g2 <- ggplot(this_df_cat_long_1) +
        geom_area(
            aes(x = start, y = value, fill = Nt)
        ) +
        scale_x_continuous("POS", labels = scales::label_number()) +
        scale_fill_manual(
            limits = nt_cat_scale_limits,
            values = nt_cat_scale_colors
        ) +
        theme_bw() +
        ggtitle(sprintf("Base distribution in %s", this_chr_name))

    g3 <- ggplot(this_df_gtf_intervals) +
        geom_bar(
            aes(x = start, y = value),
            stat = "identity"
        ) +
        facet_grid(IntervalType ~ ., scales = "free_y") +
        scale_x_continuous("POS", labels = scales::label_number()) +
        theme_bw() +
        ggtitle(sprintf("Interval distribution in %s", this_chr_name))

    g4 <- ggplot(this_df_diversity) +
        geom_area(
            aes(x = start, y = value)
        ) +
        facet_grid(DiversityIndexName ~ .) +
        scale_x_continuous("POS", labels = scales::label_number()) +
        scale_y_continuous(limits = c(0, 1)) +
        theme_bw() +
        ggtitle(sprintf("Diversity of %s", this_chr_name))

    ggsave(
        sprintf("%s/%s-nt.pdf", argv$output, this_chr_name),
        g1, width = 20, height = 5
    )
    ggsave(
        sprintf("%s/%s-nt-c.pdf", argv$output, this_chr_name),
        g2, width = 20, height = 5
    )
    ggsave(
        sprintf("%s/%s-gtf.pdf", argv$output, this_chr_name),
        g3, width = 20, height = 5
    )
    ggsave(
        sprintf("%s/%s-div.pdf", argv$output, this_chr_name),
        g4, width = 20, height = 5
    )
}

a <- parLapply(
    cl = cl,
    X = unique(df_cat$chr_name),
    plot_chr
)

parallel::stopCluster(cl)
