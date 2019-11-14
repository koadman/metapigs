# 13_1.R script                                    #
# loads unnormalized dataset                       #
# subsets, widest, aggregates                      #
# Subset to keep x number of NAs per row           #
# NAs to zeros                                     #
# runs analysis with Voom                          #

# loads normalized dataset                         #
# subsets, widest, aggregates, NAs to zeros        #
# makes long again                                 #
# plots top genes                                  #
# plots top genes with error bars                  #


library(tidyverse)
library(data.table)
library(splitstackshape)
library(limma)
library(edgeR)
library(ggplot2)
library(ggpubr)


# loads input file: tot_counts_dereplicated.csv from script 7.R (unnormalized data)
tot_counts_dereplicated <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/tot_counts_dereplicated.csv")
View(tot_counts_dereplicated)

# need to substitute cohort names that contain symbols: 
tot_counts_dereplicated[1] <- lapply(
  tot_counts_dereplicated[1], 
  gsub, 
  pattern = "Neomycin+D-scour", 
  replacement = "NeoD", 
  fixed = TRUE)
tot_counts_dereplicated[1] <- lapply(
  tot_counts_dereplicated[1], 
  gsub, 
  pattern = "Neomycin+ColiGuard", 
  replacement = "NeoC", 
  fixed = TRUE)
tot_counts_dereplicated[1] <- lapply(
  tot_counts_dereplicated[1], 
  gsub, 
  pattern = "D-scour", 
  replacement = "Dscour", 
  fixed = TRUE)

# Possible cohort names now are: 
unique(tot_counts_dereplicated$cohort)

# subset to dates: (in this case t=0 and t=2 : before and after treatmnet)
cohorts_subsets <- tot_counts_dereplicated %>% filter(
  date == "t0" | date == "t2" ,
  cohort == "Control" | cohort == "Neomycin" | cohort == "NeoD" | cohort == "NeoC" | cohort == "Dscour" | cohort == "ColiGuard" | cohort == "Mothers"
  )
NROW(cohorts_subsets)

# make widest, one row per secondary_custer + bin
cohorts_subsets_widest <- cohorts_subsets %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")

# aggregate: take average of all bins that belong to the same cluster
bbb <- cohorts_subsets_widest
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("secondary_cluster")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)
# this is how many metagenomic bins we have
NROW(ccc)


#######################################

#sort columns
ccc_sorted1 <- ccc %>% 
  select(sort(names(.)))
#secondary_cluster (metagenomic bin) column at first position
ccc_sorted2 <- ccc_sorted1 %>% 
  select("secondary_cluster", everything())


#######################################

# IMPORTANT CHANGE

# Subset to df containing rows with X number of NAs in cols

# when selecting cohorts Neomycin and Control, t0 and t2, 80 is the right number to choose 
#NCOL(ccc_sorted2)
NROW(ccc_sorted2)
#tokeep <- ccc_sorted2[rowSums(is.na(ccc_sorted2)) <= 80, ]
#NCOL(tokeep)
#NROW(tokeep)

# alternative: 
ccc_sorted2 <- data.frame(ccc_sorted2[,-1])
ccc_sorted2[is.na(ccc_sorted2)] <- 0
tokeep <- ccc_sorted2[ rowSums(cpm(ccc_sorted2)>1) >= 1, ]
NROW(tokeep)
#######################################

#replace missing values with zeros
#tokeep[is.na(tokeep)] <- 0

#######################################



# START analysis:

#######################################

# Prepare metadata

# convert countdata first column secondary_cluster (metagenomic bins) to rownames
#count_data <- data.frame(tokeep[,-1], row.names=tokeep[,1], check.names = FALSE)
count_data <- data.frame(ccc_sorted2, row.names=ccc_sorted1[,1], check.names = FALSE)


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
splitcols <- cSplit(df, "sample", "_")

#unite cols to form: sampleID (date+pig) and groupID (cohort+date)
aaa <- unite(splitcols,sampleID, 1:3, sep = "_", remove = FALSE)
#rename cols
names(aaa) <- c("sampleID", "date", "cohort", "pig")

# turn into a df (necessary to set rownames)
df <- as.data.frame(aaa)

#turn sampleID into rownames
metadata <- data.frame(df[,-1], row.names=df[,1])

#######################################

# Prepare the design matrix 

Treat <- factor(paste(metadata$cohort,metadata$date,sep="."))
# ~ 0 means that first group is NOT used an an intercept
design <- model.matrix(~0+Treat)
# rename column names of Design matrix
colnames(design) <- levels(Treat)


#######################################
#NROW(count_data)
#filtering before going into voom
#count_data <- count_data[ rowSums(cpm(count_data)>100) >= 5, ]
#NROW(count_data)


# VOOM:

count_data <- voom(count_data)

count_data_2=as.matrix(count_data)

#######################################

# Fit! 

#blocking on subject
corfit <- duplicateCorrelation(count_data_2,design,block=metadata$pig)
corfit$consensus

fit <- lmFit(count_data_2,design,block=metadata$pig,correlation=corfit$consensus)

# levels(Treat)
# "Control.t0"  "Control.t2"  "NeoC.t0"     "NeoC.t2"     "NeoD.t0"     "NeoD.t2"     "Neomycin.t0" "Neomycin.t2"

cm <- makeContrasts(
  #N1vN0 = Neomycin.t2-Neomycin.t0,
  C1vC0 = Control.t2-Control.t0,
  N0vC0 = Neomycin.t0-Control.t0,
  #H1
  N1vC1 = (Neomycin.t2-Neomycin.t0)-(Control.t2-Control.t0),
  #H1_1
  #ND1vC1 = (NeoD.t2-NeoD.t0)-(Control.t2-Control.t0),
  #H1_2
  #NC1vC1 = (NeoC.t2-NeoC.t0)-(Control.t2-Control.t0),
  #H1_3 (nothing expected)
  #N1vND1 = (Neomycin.t2-Neomycin.t0)-(NeoD.t2-NeoD.t0),
  #H1_3 (nothing expected)
  #N1vNC1 = (Neomycin.t2-Neomycin.t0)-(NeoC.t2-NeoC.t0),
  #H2
  #Ds1vC1 = (Dscour.t2-Dscour.t0)-(Control.t2-Control.t0),
  #H2_2
  #Co1vC1 = (ColiGuard.t2-ColiGuard.t0)-(Control.t2-Control.t0),
  levels = design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)


# Test hypotheses

# hypothesis 1: neomycin has an effect (neomycin group vs Control before and after treatment)
H1 <- topTable(fit2, n = 50, coef="N1vC1", adjust.method = "BH", sort.by = "logFC", resort.by = "P")
H1 <- subset(H1, adj.P.Val < 0.05)
NROW(H1)

# hypothesis 1_1: neomycin has an effect (neomycin+Dscour vs Control before and after neo treatment)
H1_1 <- topTable(fit2, n = 50, coef="ND1vC1", adjust.method = "BH", sort.by = "logFC", resort.by = "P")
H1_1 <- subset(H1_1, adj.P.Val < 0.05)
NROW(H1_1)

# hypothesis 1_2: neomycin has an effect (neomycin+ColiGuard vs Control before and after neo treatment)
H1_2 <- topTable(fit2, n = 50, coef="NC1vC1", adjust.method = "BH", sort.by = "logFC", resort.by = "P")
H1_2 <- subset(H1_2, adj.P.Val < 0.05)
NROW(H1_2)

# hypothesis 1_3: no diff should be seen between neo group and neo+Dscour before and after neo treatment
H1_3 <- topTable(fit2, n = 50, coef="N1vND1", adjust.method = "BH", sort.by = "logFC", resort.by = "P")
H1_3 <- subset(H1_3, adj.P.Val < 0.05)
NROW(H1_3)
View(H1_3)

# hypothesis 1_4: no diff should be seen between neo group and neo+ColiGuard before and after neo treatment
H1_4 <- topTable(fit2, n = 50, coef="N1vNC1", adjust.method = "BH", sort.by = "logFC", resort.by = "P")
H1_4 <- subset(H1_4, adj.P.Val < 0.05)
NROW(H1_4)
View(H1_4)



#######################################

# Plotting top genes:  
# make a plot for each of the metagenomic bins in toptable: do that using normalized data
# upload normalized data 
normalized <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/normalized.csv")

# subset
Ctrl_neo_0131_0207 <- normalized %>% filter(
  date == "t0" | date == "t2", 
  cohort == "Neomycin" | cohort == "Control" | cohort == "Neomycin+ColiGuard" | cohort == "Neomycin+D-scour"
)
View(Ctrl_neo_0131_0207)

# make widest, one row per secondary_custer + bin
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "norm_counts")

# average of all bins that belong to the same cluster
bbb <- Ctrl_neo_0131_0207_widest
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("secondary_cluster")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)

#replace missing values with zeros
ccc[is.na(ccc)] <- 0

# make long again
Ctrl_0131_0207_long <- melt(ccc, id=c("secondary_cluster"))

# (Making wide and long again was necessary to get a zero count when a pig 
# has a data point at t=0 but not for t=1 and viceversa)

#split column "variable" into date, cohort and pig. 
#split column "variable" using _ separator
split <- cSplit(Ctrl_0131_0207_long, "variable", "_")
# rename new columns
colnames(split)[colnames(split)=="variable_1"] <- "date"
colnames(split)[colnames(split)=="variable_2"] <- "cohort"
colnames(split)[colnames(split)=="variable_3"] <- "pig"

Ctrl_neo_0131_0207 <- split

###############################

# Subset to rows that match item in TOP GENES list
Ctrl_neo_0131_0207_top <- subset(Ctrl_neo_0131_0207, secondary_cluster %in% rownames(H1))

#######################################

# PLOT : Barplots with dots plus sd bars 

gg_fun <- function(parameter, dt){
  p <- ggbarplot(dt[dt$secondary_cluster == parameter, ], x = "date", y = "value", 
                 color = "cohort", palette = "Paired",
                 add = c("mean_sd", "jitter"),
                 error.plot = "upper_errorbar",
                 position = position_dodge())+
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(Ctrl_neo_0131_0207_top$secondary_cluster), gg_fun, dt = Ctrl_neo_0131_0207_top)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_top_dotsPLUSerrorbars2.pdf")
plot_list
dev.off()


####################

# subset of a cluster to play with plot: 

TTT <- Ctrl_neo_0131_0207_top %>% filter(
  secondary_cluster == "688_2"
)

# Ways to normalize y axis: 
require(scales)

p <- ggbarplot(TTT, x = "date", y = "value", 
               color = "cohort", palette = "Paired",
               add = c("mean_sd", "jitter"),
               error.plot = "upper_errorbar",
               position = position_dodge())
q <- p + scale_y_sqrt()
q


# possibilities : 
p + scale_y_sqrt()
p + scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))


#######################################


save.image(file = "~/Desktop/bins_clustering_parsing_dataframes/20191112.RData")
load("~/Desktop/bins_clustering_parsing_dataframes/20191112.RData")