
# 9.R script                                              #
# subsetting to cohort and dates of interest,             #
# make widest,                                           #
# edgeR                                                   #


############################################################################################################

# TEST HYPOTHESES using EdgeR and DESeq. you'll need to make the dataframe the widest
# start from dataframe Ctrl_neo_0131_0207 and make it the widest you can, 
# keeping unique bin + secondary_cluster as rows


############################################################################################################

# load input file (whole dataset)
total_dereplicated <- read_csv("~/Desktop/bins_clustering_parsing_dataframes/total_dereplicated.csv")
View(total_dereplicated)


# TEST HYPOTHESES: 

# subset to dates: Jan31 and Fe7
# subset to cohorts: control and neomycin
Ctrl_neo_0131_0207 <- total_dereplicated %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
)

View(Ctrl_neo_0131_0207)

#make widest, one row per secondary_custer + bin
library(tidyverse)
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")
View(Ctrl_neo_0131_0207_widest)

#join first two columns "bin + secondary_cluster" to create genome ID
Ctrl_neo_0131_0207_widest$sec_clu_bin <- paste(Ctrl_neo_0131_0207_widest$secondary_cluster, Ctrl_neo_0131_0207_widest$bin, sep="_")
#remove column bin and secondary cluster
Ctrl_neo_0131_0207_widest$bin <- NULL
Ctrl_neo_0131_0207_widest$secondary_cluster <- NULL
#bring the sec_clu_bin column to first position
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207_widest[,c(ncol(Ctrl_neo_0131_0207_widest),1:(ncol(Ctrl_neo_0131_0207_widest)-1))]

#replace NAs with zeros (yet to be determined if it's a good idea in our case)
Ctrl_neo_0131_0207_widest[is.na(Ctrl_neo_0131_0207_widest)] <- 0

# proof of why it's smart to average bins that belong to the same cluster
#z <- filter(Ctrl_neo_0131_0207, pig == "29797", date == "2017-01-31")
#View(z)
#you'll see if you sort by secondary_cluster, that some sec clusters are repeated 
#belonging to more than one bin
# it means that metabat2 failed to bin those two bins together 
#(that accordingly to dRep belong together)
#so now we need to take the avg and keep only one row for each sec cluster unique to one sample 

#I checked the above and averaging doesn t seem necessary (already averaged), but not 100% sure
# average of all bins that belong to the same cluster (works!)
#colsToAggregate <- colnames(Ctrl_neo_0131_0207_widest[,3:ncol(Ctrl_neo_0131_0207_widest)])
#aggregateBy <- c("secondary_cluster")
#dummyaggfun <- function(v, na.rm = TRUE) {
#  c(mean = mean(v, na.rm = TRUE))
#}
#Ctrl_neo_0131_0207_final <- aggregate(Ctrl_neo_0131_0207_widest[colsToAggregate], by = Ctrl_neo_0131_0207_widest[aggregateBy], FUN = dummyaggfun)


fwrite(Ctrl_neo_0131_0207_widest, file = "~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_widest.csv")








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



