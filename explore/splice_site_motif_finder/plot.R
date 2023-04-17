library(tidyverse)

s <- readr::read_csv("start.csv")
s_long <- s %>%
    tidyr::gather(key = "K", value = "V", -...1)
s_mtx <- as.data.frame(s)
rownames(s_mtx) <- s_mtx$`...1`
s_mtx$`...1` <- NULL

e <- readr::read_csv("ends.csv")
e_long <- e %>%
    tidyr::gather(key = "K", value = "V", -...1)
e_mtx <- as.data.frame(e)
rownames(e_mtx) <- e_mtx$`...1`
e_mtx$`...1` <- NULL

pheatmap::pheatmap(s_mtx[80:120,], cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap::pheatmap(e_mtx[80:120,], cluster_rows = FALSE, cluster_cols = FALSE)

ggplot(s) +
    geom_bar(
        aes(
            x = `...1`,
            y = V,
            fill = K
        ),
        stat = "identity"
    ) +
    theme_bw() +
    scale_x_continuous(
        limits = c(-20, 20)
    )
