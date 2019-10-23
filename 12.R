# 12.R script                                             #



library(tidyverse)
library(data.table)

# loads input file: tot_counts_dereplicated.csv from script 7.R
tot_counts_dereplicated <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/tot_counts_dereplicated.csv")

# subset to dates: Jan31 and Fe7
# subset to cohorts: control and neomycin
Ctrl_neo_0131_0207 <- tot_counts_dereplicated %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
)

#make widest, one row per secondary_custer + bin
library(tidyverse)
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")
View(Ctrl_neo_0131_0207_widest)

# aggregate: average of all bins that belong to the same cluster (works!)
bbb <- Ctrl_neo_0131_0207_widest
View(bbb)
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("secondary_cluster")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)

#replace missing values with zeros
ccc[is.na(ccc)] <- 0

#sort columns
ccc_sorted1 <- ccc %>% 
  select(sort(names(.)))
#secondary_cluster column at first position
ccc_sorted2 <- ccc_sorted1 %>% 
  select("secondary_cluster", everything())
View(ccc_sorted2)


# START F1000 edgeR INSTRUCTIONS:

#######################################

# Prepare metadata

# convert countdata first column secondary_cluster to rownames
mobData <- data.frame(ccc_sorted2[,-1], row.names=ccc_sorted2[,1], check.names = FALSE)

# create an empty df (purpose: fill in the metadata) based on the current data df:
empty_df = data.frame(
  sample = character())

df = rbind(
  empty_df,
  data.frame(
    sample = colnames(mobData)
  )
)

#split values in column into multiple columns by separator _
library(splitstackshape)
splitcols <- cSplit(df, "sample", "_")

#rename cols
names(splitcols) <- c("date", "cohort", "pig")
#unite cols to form: sampleID (date+pig) and groupID (cohort+date)
aaa <- unite(splitcols,sampleID, 1:3, sep = "_", remove = FALSE)
splitcols2 <- cSplit(aaa, "date", "-")
#rename cols
names(splitcols2) <- c("sampleID", "cohort", "pig", "year", "month", "day")
View(splitcols2)

# turn into a df (necessary to set rownames)
df <- as.data.frame(splitcols2)

#turn sampleID into rownames
metadata <- data.frame(df[,-1], row.names=df[,1])
View(metadata)

Group <- factor(paste(metadata$cohort, metadata$month, sep="."))
cbind(metadata, Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)

#filtering before ging into DGE list
mobData <- mobData[ rowSums(cpm(mobData)>100) >= 2, ]

y <- DGEList( counts=mobData)
y = calcNormFactors(y)
y= estimateGLMCommonDisp(y, design, verbose=TRUE)
fit <- glmQLFit(y, design)
my.contrasts <- makeContrasts(
  N1vN0 = Neomycin.2-Neomycin.1,
  C1vC0 = Control.2-Control.1,
  N0vC0 = Neomycin.1-Control.1,
  N1vC1 = (Neomycin.2-Neomycin.1)-(Control.2-Control.1),
  levels = design)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"N1vC1"])

colnames(fit)

contrast_N1vC1 <- glmLRT( fit, contrast=my.contrasts[,"N1vC1"])


topTags( contrast_N1vC1, n = 10 )

#######################################

