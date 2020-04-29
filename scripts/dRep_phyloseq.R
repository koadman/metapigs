
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
library(data.table)



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

# HEATMAP

# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 40 samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.03) > 40, TRUE)

# HEATMAP with only most abundant OTUs
# plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

sampleOrder = sort(sample_names(carbom_abund))
# method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 

# HEATMAP time - genus, family, order, etc ...
pdf("dRep_phylo_heatmap.pdf")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "gOTU_2", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster") 
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Species)")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "genus", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Genus)")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "family", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Family)")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "order", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Order)")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE,
             taxa.label = "class", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Class)")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE,
             taxa.label = "phylum", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Phylum)")
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

mdat <- mdat %>%
  dplyr::filter(!Cohort=="Sows")  %>%
  dplyr::filter(!`*collection_date`=="2017-01-31"|
                  `*collection_date`=="2017-02-01"|
                  `*collection_date`=="2017-02-03") %>%
  dplyr::select(isolation_source,PigPen)

colnames(mdat) <- c("pig","pen")
head(mdat)
mdat <- as.data.frame(mdat)
class(mdat)

mdat1 <- mdat %>%
  group_by(pig) %>%
  distinct()
NROW(mdat1)

mdat1$pen <- gsub("nan",NA,mdat1$pen)
mdat1 <- na.omit(mdat1)

# we need to keep only record of pigs that were not relocated. 
mdat2 <- setDT(mdat1)[,if(.N ==1) .SD,by=pig]


######################################################################


# SAMPLES 

# create a sample table for phyloseq 

# create a sample table for phyloseq 

sample_df <- df

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

#########################

# IMPORTANT: here is where I add the Pig Pen info: 

NROW(mdat2)
NROW(sample_df)

sample_df$pig <- as.character(sample_df$pig)
mdat2$pig <- as.character(mdat2$pig)
mdat2 <- as.data.frame(mdat2)


sample_df <- dplyr::left_join(sample_df,mdat2)
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
                                                  "t10",
                                                  "no-t"))

# reorder cohorts 
sample_df$cohort  = factor(sample_df$cohort, levels=c("Control", 
                                                      "DScour",
                                                      "ColiGuard", 
                                                      "Neomycin",
                                                      "NeoD",
                                                      "NeoC"))

rownames(sample_df) <- sample_df[,1]
# ready



# ready

######################################################################

# create phyloseq object

OTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)

###################################


# Control 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
ctrlt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
ctrlt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# DScour 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
dscourt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
dscourt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# ColiGuard 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
coligt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
coligt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# Neomycin
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
neot0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
neot8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoD 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
neoDt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
neoDt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoC 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
neoCt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
neoCt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################


ctrl <- ggarrange(ctrlt0,ctrlt8)
neo <- ggarrange(neot0,neot8)
dscour <- ggarrange(dscourt0,dscourt8)
colig <- ggarrange(coligt0,coligt8)
neoD <- ggarrange(neoDt0,neoDt8)
neoC <- ggarrange(neoCt0,neoCt8)

pdf("dRep_phylo_cohousing.pdf")
ggarrange(
  ctrl, neo, nrow=2, labels=c("A","B")
)
ggarrange(
  dscour, colig, nrow=2, labels=c("A","B")
)
ggarrange(
  neoD, neoC, nrow=2, labels=c("A","B")
)
dev.off()

########################################################################################
