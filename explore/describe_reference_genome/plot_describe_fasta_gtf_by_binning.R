library(argparser)

p <- arg_parser("")
p <- add_argument(p, "--input", help="input Parquet file")
p <- add_argument(p, "--output", help="Output Directory")
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
        gtf_intervals_strandless
    ))
df_gtf_intervals <- df %>%
    dplyr::select(
        start,
        chr_name,
        gtf_intervals_pos,
        gtf_intervals_neg,
        gtf_intervals_strandless
    )

normal_nts <- c("A", "C", "G", "T")
sm_nts <- c("a", "c", "g", "t")
other_nts <- setdiff(
    colnames(df_cat),
    c(normal_nts, sm_nts, "N", "n", "chr_name", "start")
)

parallel::clusterExport(cl, varlist = ls())

plot_chr <- function(this_chr_name){
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
        dplyr::mutate(
            NORMAL = rowSums(dplyr::select(this_df_cat, dplyr::all_of(normal_nts))),
            SOFT_MASKED = rowSums(dplyr::select(this_df_cat, dplyr::all_of(sm_nts))),
            OTHERS= rowSums(dplyr::select(this_df_cat, dplyr::all_of(other_nts)))
        ) %>%
        dplyr::select(NORMAL, SOFT_MASKED, OTHERS, N, n, start) %>%
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

    message(sprintf("Plotting %s: Plotting", this_chr_name))
    g1 <- ggplot(this_df_cat_long) +
        geom_area(
            aes(x = start, y = value, fill = Nt)
        ) +
        scale_x_continuous("POS", labels = scales::label_number()) +
        theme_bw() +
        ggtitle(sprintf("Base distribution in %s", this_chr_name))

    g2 <- ggplot(this_df_cat_long_1) +
        geom_area(
            aes(x = start, y = value, fill = Nt)
        ) +
        scale_x_continuous("POS", labels = scales::label_number()) +
        scale_color_manual() +
        theme_bw() +
        ggtitle(sprintf("Base distribution in %s", this_chr_name))

    g3 <- ggplot(this_df_gtf_intervals) +
        geom_area(
            aes(x = start, y = value)
        ) +
        facet_grid(IntervalType~., scales="free_y") +
        scale_x_continuous("POS", labels = scales::label_number()) +
        theme_bw() +
        ggtitle(sprintf("Interval distribution in %s", this_chr_name))

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
}

a <- parLapply(
    cl=cl,
    X=unique(df_cat$chr_name),
    plot_chr
)

parallel::stopCluster(cl)
