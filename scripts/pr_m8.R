library(qvalue)
m8 <- read.table('m8.txt', col.names=c('gene', 'lrt'))
m8$pvalue <- pchisq(m8$lrt, df=1, lower.tail = F)
qv <- qvalue(subset(m8, pvalue<1)$pvalue, pi0.method='bootstrap', robust=T)
m8$qvalue <- 1
m8[m8$pvalue<1, 'qvalue'] <- qv$qvalues

plot(qv)
hist(qv)

m8$detected <- m8$qvalue < 0.1
summary(m8)

write.table(m8, 'm8_results.txt')



