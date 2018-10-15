library(ape)
library(qvalue)
library(plyr)

t.r <- read.tree('concat/rooted.nwk')
plot.phylo(t.r, show.node.label=T)


## node labels = ids in the table
bs.tr <- read.tree('bs_tree.nwk')

lrt <- read.table('bs_branch_lrt.txt', col.names=c('branch', 'gene', 'lrt', 'omega2'))
lrt$pvalue <- pchisq(lrt$lrt, df=1, lower.tail = F)
qv <- qvalue(subset(lrt, pvalue<1)$pvalue, pi0.method='bootstrap', robust=T)
plot(qv)
hist(qv)
lrt$qvalue <- 1
lrt[lrt$pvalue<1,'qvalue'] <- qv$qvalues
lrt$detected <- lrt$qvalue < 0.1
write.csv(lrt, 'bs_results.txt', row.names=FALSE)

bstat <- merge(
    aggregate(detected ~ branch, data=lrt, sum),
    rename(aggregate(detected ~ branch, data=lrt, length), c('detected'='n.tests'))
)
rownames(bstat) <- bstat$branch
bstat$prop <- with(bstat, detected/n.tests)



nval.n <- bstat[bs.tr$node.label, 'n.tests']
npos.n <- bstat[bs.tr$node.label, 'detected']
nval <- c(rep(NA, length(bs.tr$tip.label)), nval.n)
npos <- c(rep(NA, length(bs.tr$tip.label)), npos.n)


prop <- npos/(nval)

ew <- nval[bs.tr$edge[,2]]
ew <- 1 + sqrt(sqrt(ew))
ew[is.na(ew) | is.infinite(ew)] <- 1



ec <- npos[bs.tr$edge[,2]]/nval[bs.tr$edge[,2]]*1000
ec[is.na(ec) | is.infinite(ec)] <- 0




col <- c('blue', rev(heat.colors(max(ec))))[ec+1]
pal <- colorRampPalette(colors=c('green','yellow', 'red'))(max(ec+1))
col <- c('blue', pal)[ec+1]

col[nval[bs.tr$edge[,2]] == 0 |
    is.na(nval[bs.tr$edge[,2]])
    ] <- 'black'

nlab <- t.r$node.label
nlab[as.numeric(nlab) >90] <-""


tip.names.df <- read.table('names.txt', stringsAsFactor=FALSE)
tip.names <- tip.names.df$V2
names(tip.names) <- tip.names.df$V1


pdf('base_tree.pdf')
plot.phylo(bs.tr,  show.tip.label=F)
nodelabels(nlab, frame='none', adj=c(1,0.5), col='brown')
dev.off()



pdf('bs_tree.pdf', width=30, height=40)
plot.phylo(bs.tr, edge.width=ew, edge.color=col)
nodelabels(paste(t.r$node.label, ifelse(is.na(nval.n), '',
                 sprintf('[%d of %d]', npos.n, nval.n)))
         , frame='none', adj=c(1,0.5))
#edgelabels()
#nodelabels(c(
#    paste(t.r$node.label, nval.n, npos.n), rep('', length(bs.tr$tip.label))
#  ), frame='none', c(-1,0))
dev.off()




## tree plot for Olga
plot(bs.tr)

bs.tr.c <- bs.tr
bs.tr.c$edge.length <- sqrt(bs.tr.c$edge.length)

el <- paste(npos[bs.tr$edge[,2]], nval[bs.tr$edge[,2]], sep='/')
el[npos[bs.tr$edge[,2]]==0] <- ""
el[is.na(npos[bs.tr$edge[,2]])] <- ""


pdf('bs_tree_vis.pdf', width=10, height=6)
plot.phylo(bs.tr.c, type='cladogram',  edge.width=ew, edge.color=col, show.tip.label=F)
edgelabels(el, frame='none')
nodelabels(nlab, frame='none', adj=c(1,0.5), col='brown')
add.scale.bar()
dev.off()

pdf('bs_tree_vis_labels.pdf', width=10, height=30)
bs.tr.c.tip.labels <- bs.tr.c
bs.tr.c.tip.labels$tip.label <- tip.names[bs.tr.c.tip.labels$tip.label]
plot.phylo(bs.tr.c.tip.labels, type='cladogram',  edge.width=ew, edge.color=col)
edgelabels(el, frame='none')
nodelabels(nlab, frame='none', adj=c(1,0.5), col='brown')
add.scale.bar()
dev.off()


pdf('bs_tree_vis_phylo.pdf', width=10, height=6)
plot.phylo(bs.tr,  edge.width=ew, edge.color=col, show.tip.label=F)
edgelabels(el, frame='none')
nodelabels(nlab, frame='none', adj=c(1,0.5), col='brown')
dev.off()

pdf('bs_tree_vis_nnum.pdf', width=20, height=30)
bs.tr.labels <- bs.tr
bs.tr.labels$tip.label <- tip.names[bs.tr.labels$tip.label]
plot.phylo(bs.tr.labels,  edge.width=ew, edge.color=col, show.tip.label=TRUE)
#edgelabels(el, frame='none')
nodelabels(bs.tr$node.label, adj=c(1,0.5))
dev.off()
