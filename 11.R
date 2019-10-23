# 11.R script                                             #
# takes non normalized data from script 8.R               #
# subsets to cohorts and dates of interest                #
# make widest                                             #
# aggregates rows: bins same cluster together             #
# sorts columns (groups)                                  #
# prepares count data and metadata                        #
# runs edgeR                                              #

# source: https://rstudio-pubs-static.s3.amazonaws.com/79395_b07ae39ce8124a5c873bd46d6075c137.html


library(tidyverse)
library(data.table)

# loads input file: tot_counts_dereplicated.csv from script 7.R
tot_counts_dereplicated <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/tot_counts_dereplicated.csv")
View(tot_counts_dereplicated)

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
View(ccc)


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
View(mobData)

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

# turn into a df (necessary to set rownames)
df <- as.data.frame(aaa)

#turn sampleID into rownames
metadata <- data.frame(df[,-1], row.names=df[,1])
View(metadata)

#######################################

# Prepare Count data

group <- paste(metadata$cohort, metadata$date, sep=".")
typeof(group)
group <- gsub('2017-01-31', '0' ,group)
group <- gsub('2017-02-07', '1' ,group)
group <- factor(group)
table(group)

# Design Matrix 
design <- model.matrix(~group)
colnames(design) <- levels(group)
design

# count data
dim(mobData)

# Filter. Atleast 100 cpm in atleast 2 samples or discard
mobData <- mobData[ rowSums(cpm(mobData)>100) >= 2, ]

# count data
dim(mobData)

# Make the DGEList
cur_DGELIST <- DGEList( counts=mobData, group=group, lib.size=colSums( mobData ) )

# Use TMM normalization
cur_DGELIST <- calcNormFactors( cur_DGELIST )

# Visualize the data
plotMDS( cur_DGELIST, method="bcv", col=as.numeric( cur_DGELIST$samples$group ), pch=20 )
legend( "bottomright", as.character( unique( cur_DGELIST$samples$group )), col=1:3, pch=20)

# Estimate dispersion
cur_DGELIST <- estimateGLMCommonDisp( cur_DGELIST, design )
cur_DGELIST <- estimateGLMTrendedDisp( cur_DGELIST, design )

cur_DGELIST <- estimateGLMTagwiseDisp( cur_DGELIST, design )
plotBCV( cur_DGELIST )

# Fitting and making comparisons
fit <- glmFit( cur_DGELIST, design )
design
  
# Make contrasts
contrast_C0_v_C1 <- glmLRT( fit, contrast=makeContrasts( Control.2017.01.31-Control.2017.02.07, levels=design ) )
contrast_N0_v_N1 <- glmLRT( fit, contrast=makeContrasts( Neomycin.2017.01.31-Neomycin.2017.02.07, levels=design ) )
contrast_C0_v_N0 <- glmLRT( fit, contrast=makeContrasts( Control.2017.01.31-Neomycin.2017.01.31, levels=design ) )

topTags(contrast_C0_v_C1)
topTags(contrast_N0_v_N1)
topTags(contrast_C0_v_N0)


# access the top genes based on a particular measure, say BH FDR < 0.05 using the decideTestsDGE function
# top genes in contrast_C0_v_C1 with a BH FDR < 0.05 and plot them with a “Smear” plot

dt_significant_contrast_C0_v_C1 <- decideTestsDGE( contrast_C0_v_C1, adjust.method="BH", p.value=0.05)
vctr_names_sig_contrast_C0_v_C1 <- rownames( cur_DGELIST )[ as.logical( dt_significant_contrast_C0_v_C1 )]
plotSmear( contrast_C0_v_C1, de.tags = vctr_names_sig_contrast_C0_v_C1 )
abline( h = c( -2, 2 ), col = "blue")

# Heatmap to get a global view of the samples
# Will highlights genes that have driven the differences

vctr_names_top <- rownames( topTags( contrast_C0_v_C1, n = 30 ) )
vctr_names_top <- c( vctr_names_top, rownames( topTags( contrast_N0_v_N1, n = 30 ) ) )
vctr_names_top <- unique( c( vctr_names_top, rownames( topTags( contrast_C0_v_N0, n = 30 ) ) ) )

vctr_sig <- as.logical( decideTestsDGE( contrast_C0_v_C1, adjust.method="BH", p.value=0.000005) )
vctr_sig <- vctr_sig | as.logical( decideTestsDGE( contrast_N0_v_N1, adjust.method="BH", p.value=0.000005) )
vctr_sig <- vctr_sig | as.logical( decideTestsDGE( contrast_C0_v_N0, adjust.method="BH", p.value=0.000005) )
vctr_names_hcl <- rownames( cur_DGELIST )[ vctr_sig ]
length( vctr_names_hcl )
head( vctr_names_hcl )

# Matrix of values
mtrx_significant <- cur_DGELIST$counts[ vctr_names_top, ]
# Colors for the groups
vctr_colors = as.factor( c( "black", "red", "green", "blue"))
vctr_colors

vctr_sample_colors <- as.character( vctr_colors[ as.numeric( cur_DGELIST$samples$group ) ] )
vctr_sample_colors

library(gplots)

# Heatmap
heatmap.2( log2( mtrx_significant + 1 ), 
           ColSideColors=vctr_sample_colors, 
           key=TRUE, 
           trace="none", 
           col=heat.colors(200), 
           scale="row" )


# More Complicated Comparisons
# Four groups

# Estimate dispersion
cur_DGELIST <- estimateGLMCommonDisp( cur_DGELIST, design )
cur_DGELIST <- estimateGLMTrendedDisp( cur_DGELIST, design )
cur_DGELIST <- estimateGLMTagwiseDisp( cur_DGELIST, design )

# Fit counts to model
fit <- glmFit( cur_DGELIST, design )

# Make contrasts comparing groups 
contrast_C0N0_v_C1N1 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Control.2017.01.31 + Neomycin.2017.01.31 )/2-( Control.2017.02.07 + Neomycin.2017.02.07 )/2, 
                                    levels=design ) )
# How did each group change over time? 
contrast_C0C1_v_N0N1 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Control.2017.01.31 + Control.2017.02.07 )/2-( Neomycin.2017.01.31 + Neomycin.2017.02.07 )/2, 
                                  levels=design ) )

# takes the difference between the successive time points, 
# e.g. it gets the deltas, and then it asks whether the deltas are different

contrast_C1N1_v_C0N0d2 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Control.2017.02.07 - Neomycin.2017.02.07 ) - ( Control.2017.01.31 - Neomycin.2017.01.31 )/2, 
                                  levels=design ) )


topTags( contrast_C0N0_v_C1N1, n = 10 )
topTags( contrast_C0C1_v_N0N1, n = 10 )
topTags( contrast_C1N1_v_C0N0d2, n = 10 )

contrast_C1N1d2_v_C0N0 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Control.2017.02.07 - Neomycin.2017.02.07 )/2 - ( Control.2017.01.31 - Neomycin.2017.01.31 ), 
                                    levels=design ) )
topTags( contrast_C1N1d2_v_C0N0, n = 10 )

contrast_C1N1_v_C0N0 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Control.2017.02.07 - Neomycin.2017.02.07 ) - ( Control.2017.01.31 - Neomycin.2017.01.31 ), 
                                    levels=design ) )
topTags( contrast_C1N1_v_C0N0, n = 10 )

contrast_N1C1_v_N0C0d2 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Neomycin.2017.02.07 - Control.2017.02.07 ) - ( Neomycin.2017.01.31 - Control.2017.01.31 )/2, 
                                  levels=design ) )
topTags( contrast_N1C1_v_N0C0d2, n = 10 )

contrast_N1C1d2_v_N0C0 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Neomycin.2017.02.07 - Control.2017.02.07 )/2 - ( Neomycin.2017.01.31 - Control.2017.01.31 ), 
                                  levels=design ) )
topTags( contrast_N1C1d2_v_N0C0, n = 10 )

contrast_N1C1_v_N0C0 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Neomycin.2017.02.07 - Control.2017.02.07 ) - ( Neomycin.2017.01.31 - Control.2017.01.31 ), 
                                    levels=design ) )
topTags( contrast_N1C1_v_N0C0, n = 10 )

contrast_C0N0_v_C1N1d2 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Control.2017.01.31 - Neomycin.2017.01.31 ) - ( Control.2017.02.07 - Neomycin.2017.02.07)/2, 
                                  levels=design ) )
topTags( contrast_C0N0_v_C1N1d2, n = 10 )

contrast_C0N0d2_v_C1N1 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Control.2017.01.31 - Neomycin.2017.01.31 )/2 - ( Control.2017.02.07 - Neomycin.2017.02.07), 
                                    levels=design ) )
topTags( contrast_C0N0d2_v_C1N1, n = 10 )

contrast_C0N0_v_C1N1 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Control.2017.01.31 - Neomycin.2017.01.31 ) - ( Control.2017.02.07 - Neomycin.2017.02.07), 
                                    levels=design ) )
topTags( contrast_C0N0_v_C1N1, n = 10 )


contrast_C1N0_v_C0N1d2 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Control.2017.02.07 - Neomycin.2017.01.31 ) - ( Control.2017.01.31 - Neomycin.2017.02.07)/2, 
                                  levels=design ) )
topTags( contrast_C1N0_v_C0N1d2, n = 10 )

contrast_C1N0d2_v_C0N1 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Control.2017.02.07 - Neomycin.2017.01.31 )/2 - ( Control.2017.01.31 - Neomycin.2017.02.07), 
                                    levels=design ) )
topTags( contrast_C1N0d2_v_C0N1, n = 10 )

contrast_C1N0_v_C0N1 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Control.2017.02.07 - Neomycin.2017.01.31 ) - ( Control.2017.01.31 - Neomycin.2017.02.07), 
                                    levels=design ) )
topTags( contrast_C1N0_v_C0N1, n = 10 )

contrast_C0C1_v_N0N1 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Control.2017.01.31 - Control.2017.02.07 ) - ( Neomycin.2017.01.31 - Neomycin.2017.02.07), 
                                  levels=design ) )
topTags( contrast_C0C1_v_N0N1, n = 10 )

contrast_C1C0_v_N1N0 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Control.2017.02.07 - Control.2017.01.31 ) - ( Neomycin.2017.02.07 - Neomycin.2017.01.31), 
                                  levels=design ) )
topTags( contrast_C1C0_v_N1N0, n = 10 )


contrast_N1N0_v_C1C0 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Neomycin.2017.02.07 - Neomycin.2017.01.31 ) - ( Control.2017.02.07 - Control.2017.01.31 ), 
                                  levels=design ) )
topTags( contrast_N1N0_v_C1C0, n = 10 )

contrast_N1N0_p_C1C0 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Neomycin.2017.02.07 - Neomycin.2017.01.31 ) + ( Control.2017.02.07 - Control.2017.01.31 ), 
                                  levels=design ) )
topTags( contrast_N1N0_p_C1C0, n = 10 )

contrast_N1pN0_p_C1pC0 <- glmLRT( fit, 
                                contrast=makeContrasts( 
                                  ( Neomycin.2017.02.07 + Neomycin.2017.01.31 ) + ( Control.2017.02.07 + Control.2017.01.31 ), 
                                  levels=design ) )
topTags( contrast_N1pN0_p_C1pC0, n = 10 )

contrast_N1pN0_v_C1pC0 <- glmLRT( fit, 
                                  contrast=makeContrasts( 
                                    ( Neomycin.2017.02.07 + Neomycin.2017.01.31 ) - ( Control.2017.02.07 + Control.2017.01.31 ), 
                                    levels=design ) )
topTags( contrast_N1pN0_v_C1pC0, n = 10 )






# Try heatmap with groups contrasts contrast_Cdiffs_v_Ndiffs and contrast_C0C1_v_N0N1
colnames(design)

# access the top genes based on a particular measure, say BH FDR < 0.05 using the decideTestsDGE function
# top genes in contrast_C0_v_C1 with a BH FDR < 0.05 and plot them with a “Smear” plot

dt_significant_contrast_C0N0_v_C1N1 <- decideTestsDGE( contrast_C0N0_v_C1N1, adjust.method="BH", p.value=0.05)
vctr_names_sig_contrast_C0N0_v_C1N1 <- rownames( cur_DGELIST )[ as.logical( dt_significant_contrast_C0N0_v_C1N1 )]
plotSmear( contrast_C0N0_v_C1N1, de.tags = vctr_names_sig_contrast_C0N0_v_C1N1 )
abline( h = c( -2, 2 ), col = "blue")

# Heatmap to get a global view of the samples
# Will highlights genes that have driven the differences

vctr_names_top <- rownames( topTags( contrast_C0N0_v_C1N1, n = 30 ) )
vctr_names_top <- c( vctr_names_top, rownames( topTags( contrast_C0C1_v_N0N1, n = 30 ) ) )

vctr_sig <- as.logical( decideTestsDGE( contrast_C0C1_v_N0N1, adjust.method="BH", p.value=0.000005) )
vctr_sig <- vctr_sig | as.logical( decideTestsDGE( contrast_C0N0_v_C1N1, adjust.method="BH", p.value=0.000005) )
vctr_names_hcl <- rownames( cur_DGELIST )[ vctr_sig ] 
length( vctr_names_hcl )
head( vctr_names_hcl )

# Matrix of values
mtrx_significant <- cur_DGELIST$counts[ vctr_names_top, ]
# Colors for the groups
vctr_colors = as.factor( c( "yellow", "red", "green", "blue"))
vctr_colors

vctr_sample_colors <- as.character( vctr_colors[ as.numeric( cur_DGELIST$samples$group ) ] )
vctr_sample_colors

library(gplots)

# Heatmap
heatmap.2( log2( mtrx_significant + 1 ), 
           ColSideColors=vctr_sample_colors, 
           key=TRUE, 
           trace="none", 
           col=heat.colors(200), 
           scale="row", 
           Colv = FALSE)

levels(group)


###############


# Try heatmap with diffs


# access the top genes based on a particular measure, say BH FDR < 0.05 using the decideTestsDGE function
# top genes in contrast_C0_v_C1 with a BH FDR < 0.05 and plot them with a “Smear” plot

dt_significant_contrast_Cdiffs_v_Ndiffs <- decideTestsDGE( contrast_Cdiffs_v_Ndiffs, adjust.method="BH", p.value=0.05)
vctr_names_sig_contrast_Cdiffs_v_Ndiffs <- rownames( cur_DGELIST )[ as.logical( dt_significant_contrast_Cdiffs_v_Ndiffs )]
plotSmear( contrast_Cdiffs_v_Ndiffs, de.tags = vctr_names_sig_contrast_Cdiffs_v_Ndiffs )
abline( h = c( -2, 2 ), col = "blue")

# Heatmap to get a global view of the samples
# Will highlights genes that have driven the differences

vctr_names_top <- rownames( topTags( contrast_Cdiffs_v_Ndiffs, n = 30 ) )
vctr_names_top <- c( vctr_names_top, rownames( topTags( contrast_Cdiffs_v_Ndiffs, n = 30 ) ) )

vctr_sig <- as.logical( decideTestsDGE( contrast_Cdiffs_v_Ndiffs, adjust.method="BH", p.value=0.000005) )
vctr_sig <- vctr_sig | as.logical( decideTestsDGE( contrast_Cdiffs_v_Ndiffs, adjust.method="BH", p.value=0.000005) )
vctr_names_hcl <- rownames( cur_DGELIST )[ vctr_sig ] 
length( vctr_names_hcl )
head( vctr_names_hcl )

# Matrix of values
mtrx_significant <- cur_DGELIST$counts[ vctr_names_top, ]
# Colors for the groups
vctr_colors = as.factor( c( "yellow", "red", "green", "blue"))
vctr_colors

vctr_sample_colors <- as.character( vctr_colors[ as.numeric( cur_DGELIST$samples$group ) ] )
vctr_sample_colors

library(gplots)

# Heatmap
heatmap.2( log2( mtrx_significant + 1 ), 
           ColSideColors=vctr_sample_colors, 
           key=TRUE, 
           trace="none", 
           col=heat.colors(200), 
           scale="row", 
           Colv = FALSE)

levels(group)
