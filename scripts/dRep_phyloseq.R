
library(readr)
library(splitstackshape)
library(tidyr)
library(phyloseq)
library(dplyr)
library(ggplot2)

library(phyloseq)
library(ggpubr)
library(pheatmap)

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(cluster)
library(circlize)
library(readxl)



setwd("~/Desktop/metapigs_dry/dRep/")
basedir = "~/Desktop/metapigs_dry/"


######################################################################

# upload bins with counts (from output of 7.R)
no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

no_reps_all$primary_cluster <- paste0(no_reps_all$secondary_cluster)
no_reps_all <- cSplit(no_reps_all,"primary_cluster","_")
no_reps_all$primary_cluster_2 <- NULL
colnames(no_reps_all)[colnames(no_reps_all)=="primary_cluster_1"] <- "primary_cluster"

colnames(no_reps_all)[colnames(no_reps_all)=="secondary_cluster"] <- "gOTU"

df <- no_reps_all
NROW(df)

######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv("/Users/12705859/Desktop/metapigs_dry/gtdbtk/gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)


######################################################################


# load metadata 
mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)
colnames(mdat)



######################################################################



# add taxa assigned by gtdbtk

df0 <- merge(df, gtdbtk_bins, by=c("pig","bin"))
NROW(df0)

######################################################################

# TAXA


taxa_mat <- df0 %>%
  group_by(gOTU) %>%
  slice(1) %>%
  dplyr::select(gOTU,phylum,class,order,family,genus,species)

taxa_mat <- cSplit(taxa_mat, "gOTU", sep="_")
taxa_mat$gOTU <- paste0(taxa_mat$gOTU_1,"_",taxa_mat$gOTU_2)
taxa_mat$gOTU_2 <- taxa_mat$gOTU
NROW(taxa_mat)
NROW(unique(taxa_mat$gOTU))

taxa_mat_df <- as.data.frame(taxa_mat)
# to matrix 
taxa_mat <- taxa_mat_df
colnames(taxa_mat)
rownames(taxa_mat) <- taxa_mat[,9]
taxa_mat[,9] <- NULL
taxa_mat <- as.matrix(taxa_mat)
# ready 

NROW(unique(rownames(taxa_mat)))
head(taxa_mat_df)

######################################################################

# gOTU  

df1 <- df0 %>%
  dplyr::select(cohort,pig,bin,date,value,gOTU)


# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"
NROW(df1)

# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,date,gOTU) %>%
  dplyr::summarise(all_bins_value = sum(value))

NROW(df2)
NROW(unique(paste0(df2$pig,df2$date)))


# assign a unique sample name 
df2$sample <- paste0(df2$date,"_",df2$pig)
# remove now pig and date (redundant)
df2$pig <- NULL
df2$date <- NULL

# long to wide 
df3 <- df2 %>%
  pivot_wider(names_from = sample, values_from = all_bins_value, values_fill = list(all_bins_value = 0)) 

# to matrix 
gOTU_mat <- as.data.frame(df3)
rownames(gOTU_mat) <- gOTU_mat[,1]
gOTU_mat[,1] <- NULL
gOTU_mat <- as.matrix(gOTU_mat)
# ready 

NROW(unique(rownames(gOTU_mat)))
NROW(unique(colnames(gOTU_mat)))


######################################################################

# SAMPLES 

# create a sample table for phyloseq 

sample_df <- df0

sample_df$sample <- paste0(sample_df$date,"_",sample_df$pig)
NROW(unique(sample_df$sample))

sample_df <- sample_df %>%
  dplyr::select(sample,pig,date,cohort) %>%
  group_by(sample) %>%
  slice(1)

sample_df$gOTU <- NULL
sample_df <- as.data.frame(sample_df)

NROW(sample_df)
head(sample_df)

unique(sample_df$cohort)

# reorder dates 
sample_df$date  = factor(sample_df$date, levels=c("t0",
                                                  "t1", 
                                                  "t2",
                                                  "t3",
                                                  "t4",
                                                  "t5",
                                                  "t6",
                                                  "t7",
                                                  "t8",
                                                  "t9",
                                                  "t10",
                                                  "tM"))

# reorder cohorts 
sample_df$cohort  = factor(sample_df$cohort, levels=c("Control", 
                                                      "DScour",
                                                      "ColiGuard", 
                                                      "Neomycin",
                                                      "NeoD",
                                                      "NeoC",
                                                      "Mothers"))

rownames(sample_df) <- sample_df[,1]
# ready


######################################################################


# create phyloseq object

OTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)

carbom <- phyloseq(OTU,TAX,samples)

############################################################################################################

# necessary to rarefy? if yes, I ll need to convert counts to integers first 
#rarecurve(t(otu_table(carbom)), step=50, cex=0.5)

# SUBSETTING phyloseq obejct

# Keep only samples to be analyzed
# carbom <- subset_samples(carbom, cohort =="ColiGuard")

carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9")))

unique(taxa_mat_df$gOTU)
# filter out bins that have not been assigned any primary cluster
carbom <- subset_taxa(carbom, (!gOTU_1 %in%  c("no")))


############################################################################################################

# NORMALIZATION 

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)

############################################################################################################

# PLOT

######################

# HEATMAP

# keep only very abundant OTUs
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.2) > 1, TRUE)

# HEATMAP with only most abundant OTUs
# plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

sampleOrder = sort(sample_names(carbom_abund))
# method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 

# HEATMAP time - genus, family, order, etc ...
pdf("dRep_phylo_heatmap.pdf")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "gOTU_1", taxa.order = "gOTU_2", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Primary clusters (dRep)")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "gOTU_2", taxa.order = "gOTU_2", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Secondary clusters (dRep)") 
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order = "gOTU_2", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Species (GTDBTK)")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "family", taxa.order = "gOTU_2", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Family (GTDBTK)") 
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "order", taxa.order = "gOTU_2", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Order (GTDBTK)") 
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "phylum", taxa.order = "gOTU_2", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Phylum (GTDBTK)") 
dev.off()



######################

# DIVERSITY 

pdf("dRep_phylo_diversity.pdf")
plot_richness(carbom, measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), 
              x="cohort", color="date")
dev.off()


######################

# ORDINATION 

carbom.ord <- ordinate(carbom, "NMDS", "bray")

pdf("dRep_phylo_ordination.pdf")
plot_ordination(carbom, carbom.ord, type="samples", color="date", #shape= "cohort", 
                title="Based on secondary clusters (99% ANI) (dRep)") + 
  geom_point(size=2) +
  facet_wrap(~cohort)
dev.off()

######################

# NETWORK ANALYSIS 


pdf("dRep_phylo_network.pdf")
ig = make_network(carbom_abund, type = "samples", distance = "bray", max.dist = 0.3)
plot_network(ig, carbom_abund, color = "date", shape = "cohort", line_weight = 0.3, 
             label = NULL, title = "sample network - Bray-Curtis distance")
dev.off()




######################################################################
######################################################################


# CO-HOUSING community structure test: 


mdat1 <- mdat %>%
  dplyr::filter(`*collection_date`=="2017-02-28") %>%
  dplyr::select(isolation_source,PigPen)


mdat1 <- as.data.frame(mdat1)
mdat1$isolation_source <- as.character(mdat1$isolation_source)

######################################################################

# SAMPLES 

# create a sample table for phyloseq 


piggies_that_stay <- mdat1$isolation_source

sample_df <- df

sample_df <- subset(sample_df, pig %in% piggies_that_stay)
NROW(sample_df)

sample_df <- as.data.frame(sample_df)
sample_df$pig <- as.character(sample_df$pig)


sample_df <- inner_join(sample_df, mdat1, by=c("pig"="isolation_source"))
NROW(sample_df)


sample_df$sample <- paste0(sample_df$date,"_",sample_df$pig)
NROW(unique(sample_df$sample))

sample_df <- sample_df %>%
  dplyr::select(sample,pig,date,cohort,PigPen) %>%
  group_by(sample) %>%
  slice(1)

sample_df$gOTU <- NULL
sample_df <- as.data.frame(sample_df)

NROW(sample_df)
head(sample_df)

#########################

# reorder dates 
sample_df$date  = factor(sample_df$date, levels=c("t0",
                                                  "t1", 
                                                  "t2",
                                                  "t3",
                                                  "t4",
                                                  "t5",
                                                  "t6",
                                                  "t7",
                                                  "t8",
                                                  "t9",
                                                  "t10"))

# reorder cohorts 
sample_df$cohort  = factor(sample_df$cohort, levels=c("Control", 
                                                      "DScour",
                                                      "ColiGuard", 
                                                      "Neomycin",
                                                      "NeoD",
                                                      "NeoC"))

rownames(sample_df) <- sample_df[,1]
# ready

######################################################################

# 
OTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)

# Control 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
carbom <- subset_samples(carbom, (date %in% c("t0")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
ctrl2 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                        title="Control t0") + 
  geom_point(size=2)
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
carbom <- subset_samples(carbom, (date %in% c("t8")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
ctrl8 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                         title="Control t8") + 
  geom_point(size=2)
######################
######################
# DScour 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
carbom <- subset_samples(carbom, (date %in% c("t0")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
DScour2 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                         title="DScour - t0") + 
  geom_point(size=2)
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
carbom <- subset_samples(carbom, (date %in% c("t8")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
DScour8 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                         title="DScour - t8") + 
  geom_point(size=2)
######################
######################
# ColiGuard 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
carbom <- subset_samples(carbom, (date %in% c("t0")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
ColiGuard2 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                           title="ColiGuard - t0") + 
  geom_point(size=2)
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
carbom <- subset_samples(carbom, (date %in% c("t8")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
ColiGuard8 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                           title="ColiGuard - t8") + 
  geom_point(size=2)
######################
######################
# Neomycin
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
carbom <- subset_samples(carbom, (date %in% c("t0")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
Neomycin2 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                              title="Neomycin - t0") + 
  geom_point(size=2)
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
carbom <- subset_samples(carbom, (date %in% c("t8")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
Neomycin8 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                              title="Neomycin - t8") + 
  geom_point(size=2)
######################
######################
# NeoD
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
carbom <- subset_samples(carbom, (date %in% c("t0")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
NeoD2 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                             title="NeoD - t0") + 
  geom_point(size=2)
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
carbom <- subset_samples(carbom, (date %in% c("t8")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
NeoD8 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                             title="NeoD - t8") + 
  geom_point(size=2)
######################
######################
# NeoC
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
carbom <- subset_samples(carbom, (date %in% c("t0")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
NeoC2 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                         title="NeoC - t0") + 
  geom_point(size=2)
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
carbom <- subset_samples(carbom, (date %in% c("t8")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
NeoC8 <- plot_ordination(carbom, carbom.ord, type="samples", color="PigPen", shape="PigPen",
                         title="NeoC - t8") + 
  geom_point(size=2)
######################
######################


ctrl <- ggarrange(
  ctrl2, ctrl8, ncol=2, common.legend = TRUE
)
DScour <- ggarrange(
  DScour2, DScour8, ncol=2, common.legend = TRUE
)
ColiGuard <- ggarrange(
  ColiGuard2, ColiGuard8, ncol=2, common.legend = TRUE
)
Neomycin <- ggarrange(
  Neomycin2, Neomycin8, ncol=2, common.legend = TRUE
)
NeoD <- ggarrange(
  NeoD2, NeoD8, ncol=2, common.legend = TRUE
)
NeoC <- ggarrange(
  NeoC2, NeoC8, ncol=2, common.legend = TRUE
)


tosave <- ggarrange(
  ctrl, Neomycin,
  DScour, ColiGuard,
  NeoD, NeoC, 
  nrow=3, 
  ncol=2,
  common.legend = FALSE,
  labels=c("A","B","C","D","E","F")
)



pdf("dRep_cohousing.pdf")
tosave
dev.off()

########################################################################################

# test distance for statistical significance: 

cohousing_df <- rbind(ctrl2$data,
      ctrl2$data,
      DScour2$data,
      DScour8$data,
      ColiGuard2$data,
      ColiGuard8$data,
      Neomycin2$data,
      Neomycin8$data,
      NeoD2$data,
      NeoD8$data,
      NeoC2$data,
      NeoC8$data)


