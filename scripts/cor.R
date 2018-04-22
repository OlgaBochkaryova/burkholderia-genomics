library(plyr)
library(ggplot2)
library(reshape2)
library(qvalue)
library(xtable)
library(broom)
library(plyr)
library(sandwich)
library(lmtest)
library(mblm)
library(car)

##source('https://bioconductor.org/biocLite.R')
##biocLite('preprocessCore')
##load package
library(preprocessCore)


m8 <- read.table('m8_all.txt', header=T)
str(m8)
m8$omega0 <- with(m8, p/(p+q))
m8.lrt <- subset(aggregate(subset(m8, select=c(OG,LRT)), by=list(m8$OG), FUN=min), select=-Group.1)

m8.par <- subset(aggregate(subset(m8, select=c(OG,omega0, omega, kappa, p0)), by=list(m8$OG), FUN=mean), select=-Group.1)
m8.agg <- merge(m8.lrt, m8.par)

m8.agg$pvalue <- pchisq(m8.agg$LRT, df=1, lower.tail = F)
qv <- qvalue(subset(m8.agg, pvalue<1)$pvalue, pi0.method='bootstrap', robust=T)
m8.agg$qvalue <- 1
m8.agg[m8.agg$pvalue<1, 'qvalue'] <- qv$qvalues

plot(qv)
hist(qv)

m8.agg$detected <- m8.agg$qvalue < 0.1
summary(m8.agg)

og.gc <- read.table('gc.txt', header=T)
og.tlen <- read.table('tlen.txt', header=T)

expr <- read.table('GEO_pseudomallei_K96423_OGs.txt', sep='\t', header=T)
expr <- expr[,c(1,3:5)]

sm1 <- read.table('common_OGs_branch2_summary.txt', header=T, stringsAsFactors=F)
sm2 <- read.table('common_OGs_branch2_summary_ad.txt', header=T)
core <- as.numeric(readLines('core_OGs.txt'))

## loosing og=19281
expr <- merge(expr, subset(sm2, select=-distance))
str(expr)

plot(expr$WT_3, expr$WT_3.1, log='xy')
cor.test(expr$WT_3, expr$WT_3.1)
cor.test(log(expr$WT_3+.00001), log(expr$WT_3.1+.00001))

expr.rpkm <-
    data.frame(sapply(expr[,2:4], function(column) 1E9 * column / expr$length / sum(column)))

## Wrong!!
## expr.rpkm <-
##     expr[,2:4] * 1E9 / expr$length / colSums(expr[,2:4])



ggplot(data=melt(subset(expr, select=c(-OG,-length))),
       aes(variable, value)) + geom_boxplot() + ylim(0, 100)

plot(expr.rpkm$WT_3, expr.rpkm$WT_3.1, log='xy')
cor.test(expr.rpkm$WT_3, expr.rpkm$WT_3.1)
cor.test(log(expr.rpkm$WT_3+.00001), log(expr.rpkm$WT_3.1+.00001))


ggplot(data=melt(expr.rpkm),
       aes(variable, value)) + geom_boxplot() + ylim(0, 100)

expr.norm <- normalize.quantiles(as.matrix(expr.rpkm))

expr.norm <- as.data.frame(expr.norm)

colnames(expr.norm) <- colnames(expr[,2:4])

ggplot(data=melt(expr.norm),
       aes(variable, value)) + geom_boxplot() + ylim(0, 10)

expr.average <- rowMeans(expr.norm)

expr.average <- cbind(expr.average, expr[,c(1,5)])

d <- merge(sm1, expr.average)
d <- merge(sm2, d)
d <- merge(d, og.gc)
d <- merge(d, og.tlen)

d$distance <- ifelse(d$distance=='changeable', NA, as.numeric(as.character(d$distance)))
d$strand.changes <- d$strand == 'changeable'
d[d$strand == 'changeable', 'strand'] <- NA
d$strand <- factor(d$strand)
d$chromosome.changes <- d$chromosome == 'changeable'
d[d$chromosome == 'changeable', 'chromosome'] <- NA
d$chromosome <- factor(d$chromosome)
d$core <- d$OG %in% core

full <- merge(d, m8.agg)

full.long <- melt(subset(full, select=c(-strand, -chromosome, -detected)),
    id.vars=c('OG'))

ggplot(
    full.long,
    aes(value)) + geom_histogram() +
    facet_wrap(~variable, scales='free_x')

range.01 <- function(v) (v-min(v, na.rm=TRUE))/(max(v, na.rm=TRUE)-min(v, na.rm=TRUE))

full.tr <- transform(full,
                  expr.average=log(expr.average + 1),
                  LRT=log(LRT+1e-6),
                  length=log(length),
                  omega0=log(omega0),
                  kappa=log(kappa),
                  GC.sd=log(GC.sd),
                  tlen=log(tlen)
                  ##distance=asin(sqrt(range.01(distance)))
                  )



full.tr.long <- melt(subset(full.tr, select=c(-strand, -chromosome, -detected)),
    id.vars=c('OG'))

ggplot(
    full.tr.long,
    aes(value)) + geom_histogram() +
    facet_wrap(~variable, scales='free_x')


full.tr$strand.leading <- full.tr$strand=='leading'
full.tr$chr1 <- full.tr$chromosome=='chr1'

high.expr.thr <- quantile(full.tr$expr.average, .8)
full.tr$high.expr <- full.tr$expr.average > high.expr.thr

na.outliers <- function(v) {
    lh <- quantile(v,probs=0.25, na.rm=T)
    uh <- quantile(v,probs=0.75, na.rm=T)
    step<- 1.5 * (uh-lh)
    out <- v < lh-step | v > lh+step
    v[out] <- NA
    v
}

my.scale <- function(df, skip=c()) {
    for (colid in seq.int(ncol(df))) {
        column <- df[[colid]]
        colname <- names(df)[colid]
        if (colname %in% skip)
            next
        if (is.numeric(column)) {
            df[colid] <- scale((column))
        }
        ##if (is.logical(column)) {
        ##    df[colid] <- scale(column)
        ##    names(df)[colid] <- paste(colname, 'TRUE', sep='.')
        ##}
        ##if (is.factor(column) & length(levels(column)) == 2) {
        ##    df[colid] <- scale(as.numeric(column))
        ##    names(df)[colid] <- paste(colname, levels(column)[2], sep='.')
        ##}
    }
    df
}


add.fl.estimate <- function(m) {
    vis.coeff(as.data.frame(summary(m)$coefficients))
}

vis.coeff <- function(d) {
    d$Estimate.str <- sprintf('%0.3f', d$Estimate)
    d$p.value <- format.pval(d$`Pr(>|t|)`, digits=3)
    d$`Std. Error` <- sprintf('%0.3f', d$`Std. Error`)
    d$`t value` <- sprintf('%0.3f', d$`t value`)
    d$p.value <- sub('e-0*(.*)', "\\\\cdot 10^{-\\1}", d$p.value)
    d$p.value <- paste('$', d$p.value, '$', sep='')
    select <- c('Estimate.str', 'Std. Error', 't value', 'p.value')
    subset(d, select=select)
}

coef2df <- function(cf) {
    cf <- rename(tidy(cf), c('estimate'='Estimate.str',
                              'p.value'='Pr(>|t|)',
                              'std.error'='Std. Error',
                              'statistic'='t value'))
    row.names(cf) <- cf$term
    vis.coeff(cf)
}


my.xtab <- function(df) print(xtable(df), type = "latex", sanitize.text.function = function(x){x})

vis.coeff <- function(d) {
    d$Estimate.str <- sprintf('%0.3f', d$Estimate)
    d$p.value <- format.pval(d$`Pr(>|t|)`, digits=3)
    d$`Std. Error` <- sprintf('%0.3f', d$`Std. Error`)
    d$`t value` <- sprintf('%0.3f', d$`t value`)
    d$p.value <- sub('e-0*(.*)', "\\\\cdot 10^{-\\1}", d$p.value)
    d$p.value <- paste('$', d$p.value, '$', sep='')
    select <- c('Estimate.str', 'Std. Error', 't value', 'p.value')
    subset(d, select=select)
}

coef2df <- function(cf) {
    cf <- rename(tidy(cf), c('estimate'='Estimate.str',
                              'p.value'='Pr(>|t|)',
                              'std.error'='Std. Error',
                              'statistic'='t value'))
    row.names(cf) <- cf$term
    vis.coeff(cf)
}

my.xtab <- function(df) print(xtable(df), type = "latex", sanitize.text.function = function(x){x})

par(mfrow=c(2,2))

full.tr.scale <- my.scale(full.tr)


## w0 ##
## core only 5 cases, removing
out.w0 <- -c(240, 466, 455, 701, 436, 998, 918, 1051, 520, 549, 1135, 939, 286)
m.w0  <- lm(omega0 ~ length + distance + strand.leading + chr1 + expr.average + GC + GC.sd + tlen, data=full.tr.scale[out.w0,])
plot(m.w0)
my.xtab(coef2df(m.w0))
summary(m.w0)

m.w0.c <- lm(omega0 ~ length + expr.average + tlen, data=full.tr.scale[out.w0,])
plot(m.w0.c)
summary(m.w0.c)
my.xtab(coef2df(m.w0.c))

## kappa ##
m.k <- lm(kappa ~ length + distance + strand.leading + chr1 + expr.average + GC + GC.sd + tlen, data=full.tr.scale)
plot(m.k)

summary(m.k)

m.k.c <- lm(kappa ~ length + chr1 + GC, full.tr.scale)
plot(m.k.c)
summary(m.k.c)

## LRT ##
m.lrt <- lm(LRT ~ length + distance + strand.leading + chr1+ expr.average + GC + GC.sd+ tlen, data=subset(full.tr.scale, LRT>0))

plot(m.lrt)

summary(m.lrt)

m.lrt.c <- lm(LRT ~ GC.sd, data=subset(full.tr.scale, LRT>0))
summary(m.lrt.c)

## detected ##
m.det <- glm(detected ~ length + distance + strand.leading + chr1 + expr.average + GC + GC.sd + core + tlen, family = binomial(link=logit), data=my.scale(full.tr, 'detected'))
residualPlots(m.det)
plot(m.det)

summary(m.det)

m.det.c <- lm(detected ~ length + tlen, data=my.scale(full.tr, 'detected'))
summary(m.det.c)

## chr1

m.chr1 <- glm(chr1 ~ length + strand.leading + expr.average + GC + GC.sd +  omega0 + tlen, data=my.scale(full.tr, 'chr1'), family=binomial(link=logit))
plot(m.chr1)
summary(m.chr1)

m.chr1.c <- lm(chr1 ~ length + expr.average + GC + tlen, data=my.scale(full.tr, 'chr1'))
plot(m.chr1.c)
summary(m.chr1.c)

## only 21 event, not significant
m.chrch <- glm(chromosome.changes ~ length + strand.leading + expr.average + GC + GC.sd + omega0 + tlen, data=my.scale(full.tr, 'crhomosome.changes'), family=binomial(link=logit))
plot(m.chrch)
summary(m.chrch)

## not significant
summary(glm(strand.leading ~ length + distance + chr1 + expr.average + GC + GC.sd  + omega0 + tlen, data=my.scale(full.tr), family=binomial(link=logit)))

## strand changes top 25 expr?
m.strch <- glm(strand.changes ~ length + distance + chr1 + expr.average + high.expr + GC + GC.sd  + omega0 + tlen, data=my.scale(full.tr), family=binomial(link=logit))
plot(m.strch)
summary(m.strch)

m.strch.c <- glm(strand.changes ~ distance + chr1, data=my.scale(full.tr), family=binomial(link=logit))
summary(m.strch.c)

m.strch.c <- lm(strand.changes ~ distance + chr1, data=my.scale(full.tr))
summary(m.strch.c)
coeftest(m.strch.c, vcov.=vcovHC)

## expr
m.expr <- lm(expr.average ~ length + chr1 + distance + strand  + GC + LRT + omega0 + tlen, data=full.tr.scale)
plot(m.expr)
summary(m.expr)
my.xtab(coef2df((m.expr)))

m.expr.c <- lm(expr.average ~ chr1  + GC  + omega0 + tlen, data=full.tr.scale)
plot(m.expr.c)
summary(m.expr.c)
my.xtab(coef2df((m.expr.c)))

m.expr.c2 <- lm(expr.average ~ chr1  + GC  + tlen, data=full.tr.scale)
plot(m.expr.c2)
summary(m.expr.c2)
my.xtab(coef2df((m.expr.c2)))



summary(lm(strand.changes ~ length + distance + expr.average + high.expr + GC + GC.sd + core + omega0, data=subset(full.tr, chr1==TRUE)))


## expression
summary(lm(expr.average ~ chr1 + distance + strand + core, data=full.tr))
hist(full.tr$expr.average)

