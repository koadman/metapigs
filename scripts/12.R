
# 12.R script                                      #
# loads unnormalized dataset                       #
# subsets, widest, aggregates, NAs to zeros        #
# runs edgeR                                       #
# loads normalized dataset                         #
# subsets, widest, aggregates, NAs to zeros        #
# makes long again                                 #
# plots top genes                                  #
# plots top genes with error bars                  #


library(tidyverse)
library(data.table)

# loads input file: tot_counts_dereplicated.csv from script 7.R (unnormalized data)
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

# aggregate: average of all bins that belong to the same cluster
bbb <- Ctrl_neo_0131_0207_widest
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




# START edgeR:

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

# turn into a df (necessary to set rownames)
df <- as.data.frame(splitcols2)

#turn sampleID into rownames
metadata <- data.frame(df[,-1], row.names=df[,1])

#######################################

# Prepare the design matrix 

Treat <- factor(paste(metadata$cohort,metadata$month,sep="."))
# ~ 0 means that first group is NOT used an an intercept
design <- model.matrix(~0+Treat)
# rename column names of Design matrix
colnames(design) <- levels(Treat)

#######################################


# Fit! 

corfit <- duplicateCorrelation(mobData,design,block=metadata$pig)
#11 warnings produced!
corfit$consensus

#filtering before going into DGE list
mobData <- mobData[ rowSums(cpm(mobData)>100) >= 2, ]

fit <- lmFit(mobData,design,block=metadata$pig,correlation=corfit$consensus)

cm <- makeContrasts(
  N1vN0 = Neomycin.2-Neomycin.1,
  C1vC0 = Control.2-Control.1,
  N0vC0 = Neomycin.1-Control.1,
  N1vC1 = (Control.2-Control.1)-(Neomycin.2-Neomycin.1),
  levels = design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
aaa <- topTable(fit2, coef="N1vC1", n=20)

View(aaa)

#######################################

# need to subtract from topTable genes that are DE in N0C0 contrast: 
bbb <- topTable(fit2, n = 50, coef="N0vC0", adjust.method = "BH")

# subtract DE genes in topTable
# based on N0vC0 topTable
aaa <- aaa[!rownames(aaa) %in% rownames(bbb), ]
NROW(aaa)
View(aaa)

# https://support.bioconductor.org/p/117545/#117557
# If you want to classify genes into groups, you could do so via decideTests() 
#with method="nestedF"; this will give you a matrix on which to categorize genes 
#based on the presence (and direction) of a response at each time point.
nested <- decideTests(fit2,method="nestedF")
eee <- as.data.frame(nested@.Data)
colnames(eee)
head(eee)


N1vN0_clusters <- eee %>% 
  rownames_to_column('secondary_cluster') %>%
  filter(N1vN0 == "1")

C1vC0_clusters <- eee %>% 
  rownames_to_column('secondary_cluster') %>%
  filter(C1vC0 == "1")

N0vC0_clusters <- eee %>% 
  rownames_to_column('secondary_cluster') %>%
  filter(N0vC0 == "1")

N1vC1_clusters <- eee %>% 
  rownames_to_column('secondary_cluster') %>%
  filter(N1vC1 == "1")



#######################################


# Estimate dispersion
#mobData= estimateGLMCommonDisp(mobData, design, verbose=TRUE)
# Calc normalization
#mobData = calcNormFactors(mobData)
#y <- DGEList( counts=mobData, group=Group)
#fit <- glmQLFit(mobData, design)
#qlf <- glmQLFTest(fit, contrast=my.contrasts[,"N1vC1"])
#contrast_N1vC1 <- glmLRT( fit, contrast=my.contrasts[,"N1vC1"])
#topTags( contrast_N1vC1, n = 20 )
#lrt <- glmLRT(fit, contrast=my.contrasts[,"N1vC1"] )
#edgeR_result <- topTags(lrt)
#deGenes <- decideTestsDGE(lrt, p=0.001)
#deGenes <- rownames(lrt)[as.logical(deGenes)]
#deGenes
#png("~/Desktop/bins_clustering_parsing_dataframes/tryplot.png", 640, 480)
#ar(mfrow=c(2,2))
#plotSmear(lrt, de.tags=deGenes)
#dev.off()

#######################################

# Plotting top genes:  
# make a plot for each of the clusters in toptable: do that using normalized data
# upload normalized data 
normalized <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/normalized.csv")

# subset
Ctrl_neo_0131_0207 <- normalized %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
)

# make widest, one row per secondary_custer + bin
library(tidyverse)
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")

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
library(splitstackshape)
split <- cSplit(Ctrl_0131_0207_long, "variable", "_")
# rename new columns
colnames(split)[colnames(split)=="variable_1"] <- "date"
colnames(split)[colnames(split)=="variable_2"] <- "cohort"
colnames(split)[colnames(split)=="variable_3"] <- "pig"

Ctrl_neo_0131_0207 <- split

###############################

# Subset to rows that match item in TOP GENES list
Ctrl_neo_0131_0207_top <- subset(Ctrl_neo_0131_0207, secondary_cluster %in% rownames(aaa))

###############################


# Perhaps I need to normalize data (like log) before plotting 


###############################

# PLOT 1 (without error bars)

library(ggplot2)
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=value, colour = cohort))+
    geom_bar(position="dodge", stat="identity") +
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(Ctrl_neo_0131_0207_top$secondary_cluster), gg_fun, dt = Ctrl_neo_0131_0207_top)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_top.pdf")
plot_list
dev.off()

###############################

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summarized
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df3 <- data_summary(Ctrl_neo_0131_0207_top, varname="value", 
                    groupnames=c("cohort", "date","secondary_cluster"))

# Convert dose to a factor variable
df3$date=as.factor(df3$date)


# PLOT 2 (with error bars)


library(ggplot2)
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=value, fill = cohort))+
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(df3$secondary_cluster), gg_fun, dt = df3)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_top_withSD.pdf")
plot_list
dev.off()

#######################################

# PLOT 3 : Barplots with data plus error bars 

library(ggpubr)
library(ggplot2)
gg_fun <- function(parameter, dt){
  p <- ggbarplot(dt[dt$secondary_cluster == parameter, ], x = "date", y = "value", 
                 add = c("mean_se", "jitter"),
                 color = "cohort", palette = c("#00AFBB", "#E7B800"),
                 position = position_dodge(0.8))+
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(Ctrl_neo_0131_0207_top$secondary_cluster), gg_fun, dt = Ctrl_neo_0131_0207_top)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_top_dataPLUSerrorbars.pdf")
plot_list
dev.off()


# Need to add lines to dots above


#######################################

save.image(file = "~/Desktop/bins_clustering_parsing_dataframes/20191029_12.RData")
load("~/Desktop/bins_clustering_parsing_dataframes/20191029_12.RData")