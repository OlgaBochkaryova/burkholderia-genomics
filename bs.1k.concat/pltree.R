#!/usr/bin/Rscript
library(ape)
library(plyr)
library(ggplot2)

told <- read.tree('../../concat/RAxML_bipartitions.bipart')
tnew <- read.tree('RAxML_bipartitions.bipart')

rt <- function(t) root(t, c('ASM95936','ASM75598'), edgelabel=TRUE)

told <- rt(told)
tnew <- rt(tnew)


tip.names.df <- read.table('../../names.txt', stringsAsFactor=FALSE)
tip.names <- tip.names.df$V2
names(tip.names) <- tip.names.df$V1

chn <- function(t) {
    t$tip.label <- tip.names[t$tip.label]
    t
}

told <- chn(told)
tnew <- chn(tnew)

plt.2 <- function(t1, t2, main1=NULL, main2=NULL, ...) {
    par(mfrow=c(1,2), mar=c(1,1,1,1))
    plot.phylo(t1, main=main1, ...)
    nodelabels(t1$node.label)
    plot.phylo(t2, main=main2, direction='leftwards', ...)
    nodelabels(t2$node.label)
}

pdf('trees.pdf', width=30, height=30)
plt.2(tnew, told, 'bootstrap, 1k', 'bootstrap, 100')
dev.off()

pdf('trees_noel.pdf', width=30, height=30)
plt.2(tnew, told, 'bootstrap, 1k', 'bootstrap, 100', use.edge.length=FALSE)
dev.off()

pdf('trees_clado_noel.pdf', width=30, height=30)
plt.2(tnew, told, 'bootstrap, 1k', 'bootstrap, 100', use.edge.length=FALSE, type='cladogram')
dev.off()

pdf('trees_clado.pdf', width=30, height=30)
plt.2(tnew, told, 'bootstrap, 1k', 'bootstrap, 100', type='cladogram')
dev.off()

bs <- function(t) na.omit(as.numeric(t$node.label))

cor.test(bs(told), bs(tnew))

theme_set(theme_bw())
pdf('bsplot.pdf')
ggplot(data.frame(old=bs(told), new=bs(tnew)), aes(old, new)) +
    geom_point(alpha=0.2) +
    labs(x='Bootstrap support value, 100 iterations', y='Bootstrap support value, 1000 iterations')
dev.off()
