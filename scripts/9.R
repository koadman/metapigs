# average of all bins that belong to the same cluster (works!)
bbb <- Ctrl_neo_0131_0207_widest
View(bbb)
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("secondary_cluster")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)
View(ccc)

fwrite(ccc, file = "~/Desktop/bins_clustering_parsing_dataframes/ccc.csv")


# edgeR

library(edgeR)
raw.data <- read.delim("~/Desktop/bins_clustering_parsing_dataframes/ccc.csv", sep=",")
View(raw.data)
d <- raw.data[,2:ncol(raw.data)]
View(d)
rownames(d) <- raw.data[,1]
View(rownames(d))

# create groups
# get the cohort IDs for each column 
invec <- colnames(d[,1:ncol(d)])
out <- rep(NA, length(invec))
for(x in c('Control', 'Neomycin')) out[grep(x, invec)] <- x
group <- out

length(group)
sum(ncol(d))

#NAs to zeros 
d[is.na(d)] <- 0

d <- DGEList(counts = d, group=group)


plotMDS.DGEList(d, main = "MDS Plot", xlim = c(-1, 1))

dim(d)
d <- calcNormFactors(d)
d
d <- estimateCommonDisp(d)
d$common.dispersion
sqrt(d$common.dispersion)
de.com <- exactTest(d)
summary(decideTestsDGE(de.com,p.value=0.01))
summary(decideTestsDGE(de.com,p.value=0.05))
com = summary(decideTestsDGE(de.com,p.value=0.05))
com_total = com[1] + com[3]

sink("~/Desktop/bins_clustering_parsing_dataframes/de_common_dispersion.txt")
topTags(de.com,n=com_total)
sink()

getPriorN(d)
d <- estimateTagwiseDisp(d, prop.used=0.5, grid.length=500)

de.tgw <- exactTest(d)
topTags(de.tgw)
summary(decideTestsDGE(de.tgw, p.value=0.01))
summary(decideTestsDGE(de.tgw, p.value=0.05))
bob = summary(decideTestsDGE(de.tgw,p.value=0.05))
total = bob[1] + bob[3]

sink("~/Desktop/bins_clustering_parsing_dataframes/de_tagwise_dispersion.txt")
topTags(de.tgw,n=total)
sink()


# DESeq

library("DESeq")
countsTable <- read.delim("~/Desktop/bins_clustering_parsing_dataframes/ccc.csv", header=TRUE, stringsAsFactors=TRUE, sep=",")
rownames(countsTable) <- countsTable$secondary_cluster
countsTable <- countsTable[,-1]
View(countsTable)

# create groups
# get the cohort IDs for each column 
invec <- colnames(d[,1:ncol(d)])
out <- rep(NA, length(invec))
for(x in c('Control', 'Neomycin')) out[grep(x, invec)] <- x
conds <- out

#I get stopped here because my data hasn t got integers
cds <- newCountDataSet (countsTable, conds)



# baySeq

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("baySeq")
library(baySeq)

data <- read.delim("~/Desktop/bins_clustering_parsing_dataframes/ccc.csv", sep=",", header = TRUE)
View(data)
rownames(data) <- data$secondary_cluster
countsTable <- countsTable[,-1]
View(countsTable)

#groups
invec <- colnames(data[,1:ncol(data)])
out <- rep(NA, length(invec))
for(x in c('Control', 'Neomycin')) out[grep(x, invec)] <- x
replicates <- out

#bayseq needs replicates so skip bayseq



# let's try BEEM!
# https://github.com/lch14forever/BEEM
devtools::install_github('csb5/beem')



