library('reshape2')
library('dplyr')
stimuli <- c(
    'bat',
    'rooster',
    'crab',
    'dragon',
    'ladybug',
    'tortoise',
    'whale',
    'ant',
    'butterfly',
    'bear',
    'bee',
    'camel',
    'cat',
    'caterpillar',
    'cow',
    'deer',
    'dog',
    'duck',
    'eagle',
    'elephant',
    'fish',
    'fly',
    'fox',
    'frog',
    'giraffe',
    'gorilla',
    'goat',
    'horse',
    'kangaroo',
    'lion',
    'lobster',
    'monkey',
    'mouse',
    'ostrich',
    'owl',
    'peacock',
    'penguin',
    'pig',
    'rabbit',
    'racoon',
    'rhinoceros',
    'seahorse',
    'sheep',
    'snail',
    'snake',
    'spider',
    'squirrel',
    'swan',
    'tiger',
    'zebra',
    'cymbals',
    'shell',
    'ski',
    'slide',
    'tractor',
    'web',
    'yo-yo',
    'airplane',
    'anchor',
    'barrel',
    'basket',
    'bell',
    'cake',
    'cannon',
    'chain',
    'church',
    'cigar',
    'clock',
    'crown',
    'dress',
    'drum',
    'fence',
    'flute',
    'glove',
    'guitar',
    'gun',
    'hammer',
    'harp',
    'helicopter',
    'iron',
    'kite',
    'ladder',
    'mitten',
    'motorbike',
    'nut',
    'piano',
    'pram',
    'skirt',
    'sledge',
    'snowman',
    'suitcase',
    'swing',
    'toaster',
    'train',
    'trumpet',
    'vase',
    'violin',
    'wheel',
    'whistle',
    'windmill'
)

d.wide <- read.csv('allfeatures_allwords.csv', na.strings = '#N/A', stringsAsFactors = FALSE)
names(d.wide)[names(d.wide) == 'yo.yo'] <- 'yo-yo'
names(d.wide)[names(d.wide) == 'raccoon'] <- 'racoon'
names(d.wide)[names(d.wide) == 'aeroplane'] <- 'airplane'
names(d.wide)[names(d.wide) == 'ladybird'] <- 'ladybug'
names(d.wide)[names(d.wide) == 'motorcycle'] <- 'motorbike'
names(d.wide)[names(d.wide) == 'sled'] <- 'sledge'
s <- names(d.wide)[3:ncol(d.wide)]
z <- stimuli %in% s
stimuli[!z]
z <- s %in% stimuli
z <- c(T,T,z)
d.wide <- d.wide[,z]
d <- reshape2::melt(d.wide, c('type','feature'), variable.name = 'concept', value.name = 'hasfeature')
head(d)
d <- d %>%
    group_by(feature) %>%
    mutate(hasfeature = (hasfeature - mean(hasfeature))) %>%
    ungroup()

d$concept <- as.factor(d$concept)
d$type <- as.factor(d$type)
d$feature <- as.factor(d$feature)
summary(d)
str(d)

M <- matrix(0, nrow = 409, ncol = 100)
for (i in 1:length(stimuli)) {
    s <- stimuli[i]
    M[,i] <- filter(d, concept == s)$hasfeature
}
z <- !apply(M,1,function(x) all(is.nan(x)))
M <- M[z,]
colnames(M) <- stimuli
# Compute cosine similarity
# cosine() calculates a similarity matrix between all **column vectors** of a matrix.
D.cos <- as.dist(lsa::cosine(M))
write.csv(file = "../cosine.csv", x=lsa::cosine(M))
#write.csv(file = "C:/Users/mbmhscc4/MATLAB/ECOG/naming/data/targets/similarity/semantic/Dilkinea/labels.csv", x=stimuli)
sqrt_truncate_r <- function(S, tau) {
# @S: n x n Similarity matrix
# @r: rank tuning parameter
# info: finds square root of S using eigen decompostion and truncates to
# rank r
    s <- svd(S);
    n <- length(s$d)
    for (r in 1:n) {
        C = s$u[,1:r] %*% diag(sqrt(s$d[1:r]));
        Sz <- C %*% t(C)
        objfunc = norm(S-Sz, type='F') / norm(S, type='F');
        if (objfunc <= tau) {
            return(list(C=C,Sz=Sz,r=r))
        }
    }
}

D <- dist(t(M),method = 'manhattan')
C <- cmdscale(D,k=3)

png('lowrank2.png', res = 300, width = 3580, height = 1500)
plot(h.lowrank)
dev.off()
png('fullrank2.png', res = 300, width = 3580, height = 1500)
plot(h.fullrank)
dev.off()

u <- svd(x = t(M))
tmp <- (u$u %*% diag(u$d))[,1:3]
dd$svd.x <- tmp[,1]
dd$svd.y <- tmp[,2]
dd$svd.z <- tmp[,3]

d.next <- read.csv('C:/Users/mbmhscc4/MATLAB/ECOG/naming/data/targets/embedding/semantic/NEXT/NEXT_CK_KIND_5D.csv', na.strings = '#N/A', stringsAsFactors = FALSE, header=F)
dd <- data.frame(
    x = C[,1],
    y = C[,2],
    z = C[,3],
    svd.x = tmp[,1],
    svd.y = tmp[,2],
    svd.z = tmp[,3],
    next.x = d.next[,1],
    next.y = d.next[,2],
    next.z = d.next[,3],
    label = stimuli)

p <- plotly::plot_ly(dd, x = ~next.x, y = ~next.y, z = ~next.z, text = ~label) %>%
    plotly::add_markers() %>%
    plotly::layout(scene = list(xaxis = list(title = 'X'),
                        yaxis = list(title = 'Y'),
                        zaxis = list(title = 'Z')))

print(p)


d.next <- read.csv('C:/Users/mbmhscc4/MATLAB/ECOG/naming/data/targets/embedding/semantic/NEXT/NEXT_CK_KIND_5D.csv', na.strings = '#N/A', stringsAsFactors = FALSE)

D.leuven <- as.matrix(read.csv("C:/Users/mbmhscc4/MATLAB/ECOG/naming/data/targets/similarity/semantic/Leuven/cosine.csv",header=F))
L.leuven <- read.csv("C:/Users/mbmhscc4/MATLAB/ECOG/naming/data/targets/similarity/semantic/Leuven/labels.txt",header=F)[[1]]
colnames(D.leuven) <- L.leuven
hc.leuven <- hclust(as.dist(1-D.leuven))
png(filename = 'Leuven_fullrank_hclust.png', width=500, height=1200)
plot(as.phylo(hc.leuven), cex = 0.9, label.offset = 0.01, main='Leuven Full-rank Hierarchical Clustering\nof Cosine Distance')
dev.off()

hc.dilkina <- hclust(as.dist(1-D.cos))
png(filename = 'Dilkina_fullrank_hclust.png', width=500, height=1200)
plot(as.phylo(hc.dilkina), cex = 0.9, label.offset = 0.01, main='Dilkina Full-rank Hierarchical Clustering\nof Cosine Distance')
dev.off()
