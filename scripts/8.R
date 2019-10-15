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

library("DESeq")
d <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/ccc.csv", row.names=1)
dim(d)
View(d)

# first col to rownames
#rownames(df) <- df[,1]
#df[,1] <- NULL
#View(df)

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

# create the DGEList object and calculate variance
ddd <- DGEList(counts = d, group=group)


dim(ddd)

ddd <- calcNormFactors(ddd)
ddd$counts
ddd$samples

plotMDS(ddd, main = "MDS Plot for Li Data", xlim = c(-1, 1))

ddd <- calcNormFactors(ddd)
ddd <- estimateCommonDisp(ddd,verbose=TRUE)

ddd <- estimateTagwiseDisp(ddd)
de.tgw <- exactTest(ddd)



de.tgw <- exactTest (ddd,common.disp=FALSE)
topTags(de.tgw)


summary(decideTestsDGE(de.tgw,p.value=0.05))
View(decideTestsDGE(de.tgw, p.value=0.05))



