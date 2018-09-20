library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(scales)

d <- read.table('bs.txt', sep=' ', header=TRUE) %>%
    rename(bs100=binning.fulltree, bs1k=old_tree_new_bs)

pdf('bootstrap_100_vs_1k.pdf', width=6, height=5)
ggplot(d, aes(bs100, bs1k)) +
    geom_bin2d(bins=50) +
    labs(x='Bootstrap support values, 100 iterations', y='Bootstrap support values, 1000 iterations') +
    scale_fill_distiller(limits=c(0, 1500), trans='sqrt', oob=squish, palette='Spectral')
dev.off()

with(d, cor.test(bs100, bs1k))
with(d, table(bs100 >= 95, bs1k >= 95, dnn=c('bs100', 'bs1k')))
