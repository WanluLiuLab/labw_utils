library(tidyverse)


csv_filename <- "/home/yuzj/Documents/cpptetgs_experimental/test_data/fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz.summary.d/summary.csv.xz"

df <- readr::read_csv(csv_filename) %>%
    dplyr::mutate(
        normal = A + G + C + T,
        SM = a + g + c + t
    ) %>%
    dplyr::select(!c(A, G, C, T, a, g, c, t))

df_cat <- df %>%
    dplyr::select(normal, SM, N, n, start, chr_name) %>%
    tidyr::gather(
        key = "Nt",
        value = "value",
        -start,
        -chr_name
    ) %>%
    dplyr::mutate(Nt = factor(Nt,levels=c("n", "N", "SM", "normal")))

df_gtf_intervals <- df %>%
    dplyr::group_by(chr_name) %>%
    dplyr::mutate(
        gtf_intervals=gtf_intervals/max(gtf_intervals)*100000
    ) %>%
    dplyr::ungroup()

# "normal", "SM", "N", "n"
p <- ggplot(df_cat) +
    geom_area(aes(x=start, y=value, fill=Nt)) +
    geom_point(
        data=df_gtf_intervals,
        aes(x=start, y=gtf_intervals),
        size=0.2
    ) +
    theme_bw() +
    facet_grid(chr_name~.)

ggsave(
    "1.png",
    p,
    width=10,
    height=50,
    limitsize = FALSE
)
