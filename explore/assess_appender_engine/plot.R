library(tidyverse)

df <- readr::read_tsv(
    "bench_result.tsv",
    show_col_types = FALSE,
    quote = "'"
) %>%
    dplyr::mutate(THROUGHPUT = 1 / TIME_SPENT * 1000)

p <- ggplot(df, aes(x = THREAD_NUM, y = THROUGHPUT)) +
    stat_summary(aes(color = APPENDER_CLASS_NAME), fun = "mean", geom = "line") +
    scale_y_continuous(trans = "log10", n.breaks = 20, labels = scales::label_number(scale_cut = scales::cut_si(unit = "l"))) +
    facet_grid(~BUFF_SIZE) +
    theme_bw()
ggsave("b.png", p, width = 12, height = 8)
