

# 10.R script                                             #
# sort columns                                            #
# edgeR from new manual + grouping working                #


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


# EdgeR from new manual:

# convert countdata first column secondary_cluster to rownames
count_data <- data.frame(Ctrl_neo_0131_0207_widest_sorted2[,-1], row.names=Ctrl_neo_0131_0207_widest_sorted2[,1], check.names = FALSE)
View(count_data)

# Extract metadata

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

# START F1000 EDGER INSTRUCTIONS
#CellType is date 
#Status is cohort
group <- paste(metadata$date, metadata$cohort, sep=".")
group <- factor(group)
table(group)

# count data
colnames(count_data) 
rownames(count_data)
dim(count_data)
head(count_data)

# is it normal that my lib size is one? 
library(edgeR)
y <- DGEList(count_data[,-1], group=group,
             genes=count_data[,1,drop=FALSE])
#?? options(digits=3)
calcNormFactors(y)

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
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

y <- estimateDisp(y, design, robust=TRUE)


fit <- glmQLFit(y, design, robust=TRUE)
# can 't make dispersion values! it says: 
# Error in glmQLFit.DGEList(y, design, robust = TRUE) : 
# No dispersion values found in DGEList object.

