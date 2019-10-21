

# 10.R script                                             #
# sorts columns                                           #
# edgeR from new manual + grouping working                #
# dispersion estimated with min.row.sum=0, then it works  #


# load input file (normalized and subset dataset from 8.R)
Ctrl_neo_0131_0207_widest <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_widest.csv", 
                                      check.names = FALSE)
View(Ctrl_neo_0131_0207_widest)

#sort columns
Ctrl_neo_0131_0207_widest_sorted1 <- Ctrl_neo_0131_0207_widest %>% 
  select(sort(names(.)))
#secondary_cluster column at first position
Ctrl_neo_0131_0207_widest_sorted2 <- Ctrl_neo_0131_0207_widest_sorted1 %>% 
  select("secondary_cluster", everything())
View(Ctrl_neo_0131_0207_widest_sorted2)


# START F1000 edgeR INSTRUCTIONS:

#######################################

# Prepare metadata

# convert countdata first column secondary_cluster to rownames
count_data <- data.frame(Ctrl_neo_0131_0207_widest_sorted2[,-1], row.names=Ctrl_neo_0131_0207_widest_sorted2[,1], check.names = FALSE)
View(count_data)

# create an empty df (purpose: fill in the metadata) based on the current data df:
empty_df = data.frame(
  sample = character())

df = rbind(
  empty_df,
  data.frame(
    sample = colnames(count_data)
  )
)

#split values in column into multiple columns by separator _
library(splitstackshape)
splitcols <- cSplit(df, "sample", "_")

#rename cols
names(splitcols) <- c("date", "cohort", "pig")
#unite cols to form: sampleID (date+pig) and groupID (cohort+date)
aaa <- unite(splitcols,sampleID, 1:3, sep = "_", remove = FALSE)

# turn into a df (necessary to set rownames)
df <- as.data.frame(aaa)

#turn sampleID into rownames
metadata <- data.frame(df[,-1], row.names=df[,1])
View(metadata)
fwrite(metadata, file = "~/Desktop/bins_clustering_parsing_dataframes/metadata.csv")

#######################################

# Prepare Count data

group <- paste(metadata$cohort, metadata$date, sep=".")
typeof(group)
group <- gsub('-', '.' ,group)
group <- factor(group)
table(group)

# count data
colnames(count_data) 
rownames(count_data)
dim(count_data)

library(edgeR)
y <- DGEList(count_data[,], group=group,
             genes=count_data[,1,drop=FALSE])

y$samples

# Filtering to remove low counts
keep <- rowSums(cpm(y) > 0.5) >= 2
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalization for composition bias (what does it do?)
y <- calcNormFactors(y)
y$samples

levels(group)
colors <- rep(c("darkgreen", "red", "blue", "yellow"), 2)
plotMDS(y, col=colors[group], pch=1)
legend("topright", legend=levels(group), col=colors, ncol=2)
plotMDS.DGEList(y, main = "MDS Plot", col=colors[group], pch=1)

# Design Matrix 
design <- model.matrix(~group)
colnames(design) <- levels(group)
design

library(statmod)
y <- estimateDisp(y, design, robust=TRUE, min.row.sum=0)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)
summary(fit$df.prior)

table(group)
# Control.2017.01.31  Control.2017.02.07 Neomycin.2017.01.31 Neomycin.2017.02.07 

B.LvsP_1 <- makeContrasts(Control.2017.01.31-Control.2017.02.07, levels=design)
B.LvsP_2 <- makeContrasts(Neomycin.2017.01.31-Neomycin.2017.02.07, levels=design)
B.LvsP_3 <- makeContrasts(Control.2017.01.31-Neomycin.2017.01.31, levels=design)
B.LvsP_4 <- makeContrasts(Control.2017.02.07-Neomycin.2017.02.07, levels=design)

res_1 <- glmQLFTest(fit, contrast=B.LvsP_1)
res_2 <- glmQLFTest(fit, contrast=B.LvsP_2)
res_3 <- glmQLFTest(fit, contrast=B.LvsP_3)
res_4 <- glmQLFTest(fit, contrast=B.LvsP_4)

topTags(res_1)
topTags(res_2)
topTags(res_3)
topTags(res_4)

