library(tidyverse)

df <- readr::read_csv("out.csv")

ggplot(df) + geom_boxplot(aes(x=FN, y=TIME)) + facet_grid(POW~., scales="free")
