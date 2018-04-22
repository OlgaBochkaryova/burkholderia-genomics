library(ape)
## root thre tree
t <- read.tree('concat/RAxML_bipartitions.bipart')

plot.phylo(t, use.edge.length=F, show.node.label=T)


plot.phylo(t,show.node.label=T)

t.r <- root(t, edgelabel=T, outgroup=c('ASM75598', 'ASM95936'))
write.tree(t.r, 'concat/rooted.nwk')
plot.phylo(t.r, show.node.label=T)
##
