library('dplyr')
library('ggplot2')

S <- as.matrix(read.csv(file = 'data/targets/similarity/semantic/Dilkinea/cosine.csv', header = FALSE))
labs <- read.csv(file = 'data/targets/similarity/semantic/Dilkinea/labels.txt', header = FALSE, stringsAsFactors = FALSE)[[1]]
rownames(S) <- labs
colnames(S) <- labs

S.ani <- S[1:50,1:50]
S.ina <- S[51:100,51:100]
D.ani <- as.dist(1-S.ani)
D.ina <- as.dist(1-S.ina)

hc.ani <- hclust(D.ani)
hc.ina <- hclust(D.ina)

plot(hc.ani)
plot(hc.ina)

Splits <- data.frame(
    index=1:100,
    label=labs,
    type=factor(rep(c(1,2),c(50,50)),levels=1:2,labels=c('animate','inanimate')),
    split=numeric(100),
    step=numeric(100)
)
da <- D.ani
di <- D.ina
for (i in 1:25) {
    x <- hclust(da)
    y <- x$labels[-x$merge[1,]]
    Splits$split[Splits$label==y[1]] <- 1
    Splits$split[Splits$label==y[2]] <- 2
    Splits$step[Splits$label==y[1]] <- i
    Splits$step[Splits$label==y[2]] <- i
    ma <- as.matrix(da)
    da <- as.dist(ma[x$merge[1,],x$merge[1,]])

    x <- hclust(di)
    y <- x$labels[-x$merge[1,]]
    Splits$split[Splits$label==y[1]] <- 1
    Splits$split[Splits$label==y[2]] <- 2
    Splits$step[Splits$label==y[1]] <- i
    Splits$step[Splits$label==y[2]] <- i
    mi <- as.matrix(di)
    di <- as.dist(mi[x$merge[1,],x$merge[1,]])
}

tabulate(Splits$split)

tmp <- Splits %>%
    filter(type=='animate',split==1) %>%
    arrange(step)
S.aniS1 <- S[tmp$index,tmp$index]
tmp <- Splits %>%
    filter(type=='animate',split==2) %>%
    arrange(step)
S.aniS2 <- S[tmp$index,tmp$index]
z <- lower.tri(x = matrix(0,ncol=25,nrow=25))
cat(sprintf('Animate Split Correlation\n'))
cat(sprintf('-------------------------\n'))
cat(sprintf('Pearson: %.4f\n', cor(S.aniS1[z],S.aniS2[z],method = 'pearson')))
cat(sprintf('Spearman: %.4f\n', cor(S.aniS1[z],S.aniS2[z],method = 'spearman')))
cat(sprintf('Kendall: %.4f\n', cor(S.aniS1[z],S.aniS2[z],method = 'kendall')))
cat(sprintf('\n'))

tmp <- Splits %>%
    filter(type=='inanimate',split==1) %>%
    arrange(step)
S.inaS1 <- S[tmp$index,tmp$index]
tmp <- Splits %>%
    filter(type=='inanimate',split==2) %>%
    arrange(step)
S.inaS2 <- S[tmp$index,tmp$index]
z <- lower.tri(x = matrix(0,ncol=25,nrow=25))
cat(sprintf('Inanimate Split Correlation\n'))
cat(sprintf('---------------------------\n'))
cat(sprintf('Pearson: %.4f\n', cor(S.inaS1[z],S.inaS2[z],method = 'pearson')))
cat(sprintf('Spearman: %.4f\n', cor(S.inaS1[z],S.inaS2[z],method = 'spearman')))
cat(sprintf('Kendall: %.4f\n', cor(S.inaS1[z],S.inaS2[z],method = 'kendall')))
cat(sprintf('\n'))

tmp <- Splits %>%
    filter(split==1) %>%
    arrange(step,type)
S.S1 <- S[tmp$index,tmp$index]
tmp <- Splits %>%
    filter(split==2) %>%
    arrange(step,type)
S.S2 <- S[tmp$index,tmp$index]
z <- lower.tri(x = matrix(0,ncol=50,nrow=50))
cat(sprintf('Overall Split Correlation\n'))
cat(sprintf('-------------------------\n'))
cat(sprintf('Pearson: %.4f\n', cor(S.S1[z],S.S2[z],method = 'pearson')))
cat(sprintf('Spearman: %.4f\n', cor(S.S1[z],S.S2[z],method = 'spearman')))
cat(sprintf('Kendall: %.4f\n', cor(S.S1[z],S.S2[z],method = 'kendall')))
cat(sprintf('\n'))

# write(as.character(Splits$label), file = "filters/DilkineaSplit/labels.txt", ncolumns = 1)
# write(as.numeric(Splits$split-1),file = "filters/DilkineaSplit/filter.txt",ncolumns = 1)

# # Random splits
# x <- rep(c(1,2),25)
# niter <- 1000
# SplitsRandomCor <- data.frame(
#     iter=rep(1:niter,3),
#     type=factor(rep(c(1,2,3),c(niter,niter,niter)),levels=1:3,labels=c('overall','animate','inanimate')),
#     pearson=numeric(niter*3),
#     spearman=numeric(niter*3),
#     kendall=numeric(niter*3)
# )
# for (i in 1:niter) {
#     SplitsRandom <- Splits
#     SplitsRandom$split <- c(sample(x,50),sample(x,50))
#
#     tmp <- SplitsRandom %>%
#         filter(split==1) %>%
#         arrange(step,type)
#     S.S1 <- S[tmp$index,tmp$index]
#     tmp <- SplitsRandom %>%
#         filter(split==2) %>%
#         arrange(step,type)
#     S.S2 <- S[tmp$index,tmp$index]
#     z <- lower.tri(x = matrix(0,ncol=50,nrow=50))
#     zi <- SplitsRandomCor$iter == i & SplitsRandomCor$type == 'overall'
#     SplitsRandomCor$pearson[zi] <- cor(S.S1[z],S.S2[z],method = 'pearson')
#     SplitsRandomCor$spearman[zi] <- cor(S.S1[z],S.S2[z],method = 'spearman')
#     SplitsRandomCor$kendall[zi] <- cor(S.S1[z],S.S2[z],method = 'kendall')
#
#     tmp <- SplitsRandom %>%
#         filter(type=='animate',split==1) %>%
#         arrange(step)
#     S.aniS1 <- S[tmp$index,tmp$index]
#     tmp <- SplitsRandom %>%
#         filter(type=='animate',split==2) %>%
#         arrange(step)
#     S.aniS2 <- S[tmp$index,tmp$index]
#     z <- lower.tri(x = matrix(0,ncol=25,nrow=25))
#     zi <- SplitsRandomCor$iter == i & SplitsRandomCor$type == 'animate'
#     SplitsRandomCor$pearson[zi] <- cor(S.aniS1[z],S.aniS2[z],method = 'pearson')
#     SplitsRandomCor$spearman[zi] <- cor(S.aniS1[z],S.aniS2[z],method = 'spearman')
#     SplitsRandomCor$kendall[zi] <- cor(S.aniS1[z],S.aniS2[z],method = 'kendall')
#
#     tmp <- SplitsRandom %>%
#         filter(type=='inanimate',split==1) %>%
#         arrange(step)
#     S.inaS1 <- S[tmp$index,tmp$index]
#     tmp <- SplitsRandom %>%
#         filter(type=='inanimate',split==2) %>%
#         arrange(step)
#     S.inaS2 <- S[tmp$index,tmp$index]
#     z <- lower.tri(x = matrix(0,ncol=25,nrow=25))
#     zi <- SplitsRandomCor$iter == i & SplitsRandomCor$type == 'inanimate'
#     SplitsRandomCor$pearson[zi] <- cor(S.inaS1[z],S.inaS2[z],method = 'pearson')
#     SplitsRandomCor$spearman[zi] <- cor(S.inaS1[z],S.inaS2[z],method = 'spearman')
#     SplitsRandomCor$kendall[zi] <- cor(S.inaS1[z],S.inaS2[z],method = 'kendall')
# }
# cat(sprintf('Animate Split Correlation (max of random splits)\n'))
# cat(sprintf('-------------------------\n'))
# ix <- which.max(filter(SplitsRandomCor,type=='animate')$pearson)
# cat(sprintf('Pearson: %.4f\n', filter(SplitsRandomCor,type=='animate')[ix,'pearson']))
# cat(sprintf('Spearman: %.4f\n', filter(SplitsRandomCor,type=='animate')[ix,'spearman']))
# cat(sprintf('Kendall: %.4f\n', filter(SplitsRandomCor,type=='animate')[ix,'kendall']))
# cat(sprintf('\n'))
#
# cat(sprintf('Inanimate Split Correlation (max of random splits)\n'))
# cat(sprintf('---------------------------\n'))
# ix <- which.max(filter(SplitsRandomCor,type=='inanimate')$pearson)
# cat(sprintf('Pearson: %.4f\n', filter(SplitsRandomCor,type=='inanimate')[ix,'pearson']))
# cat(sprintf('Spearman: %.4f\n', filter(SplitsRandomCor,type=='inanimate')[ix,'spearman']))
# cat(sprintf('Kendall: %.4f\n', filter(SplitsRandomCor,type=='inanimate')[ix,'kendall']))
# cat(sprintf('\n'))
#
# cat(sprintf('Overall Split Correlation (max of random splits)\n'))
# cat(sprintf('-------------------------\n'))
# ix <- which.max(filter(SplitsRandomCor,type=='overall')$pearson)
# cat(sprintf('Pearson: %.4f\n', filter(SplitsRandomCor,type=='overall')[ix,'pearson']))
# cat(sprintf('Spearman: %.4f\n', filter(SplitsRandomCor,type=='overall')[ix,'spearman']))
# cat(sprintf('Kendall: %.4f\n', filter(SplitsRandomCor,type=='overall')[ix,'kendall']))
# cat(sprintf('\n'))
#
