
library(readr)
library(splitstackshape)
library(tidyr)
library(phyloseq)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/metapigs_dry/gtdbtk")
basedir = "~/Desktop/metapigs_dry/"


# input files: 
# gtdbtk_bins_completeTaxa
# no_reps_all.csv (BINS COUNTS)
# Metagenome.environmental_20190308_2.xlsx (metadata, necessary for last part, co-housing)


# OUTPUTS:

# gt_phylo_barplot.pdf
# gt_phylo_heatmap.pdf
# gt_phylo_ordination.pdf
# gt_phylo_diversity.pdf
# gt_phylo_network.pdf
# gt_cohousing.pdf


######################################################################

# counts data 

no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)


######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv("gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################

# merge info 

NROW(gtdbtk_bins)
NROW(no_reps_all)
head(gtdbtk_bins)
head(no_reps_all)
df <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df)[colnames(df) == 'node'] <- 'gOTU'
df$gOTU <- paste0("gOTU_",df$gOTU)

NROW(unique(df$gOTU))
NROW(df)


######################################################################

# TAXA


taxa_mat <- df %>%
  dplyr::select(gOTU,species,genus,family,order,class,phylum,domain) %>%
  group_by(gOTU) %>%
  slice(1)

NROW(taxa_mat)
NROW(unique(taxa_mat$gOTU))

taxa_mat_df <- as.data.frame(taxa_mat)
# to matrix 
taxa_mat <- taxa_mat_df
rownames(taxa_mat) <- taxa_mat[,1]
taxa_mat[,1] <- NULL
taxa_mat <- as.matrix(taxa_mat)
# ready 

NROW(unique(rownames(taxa_mat)))
head(taxa_mat_df)


######################################################################

# gOTU  

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","gOTU")
df1 <- df[ , (names(df) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

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

#subset to what you want & remove phyla == NA
unique_phyla <- unique(taxa_mat_df$phylum)
unique_phyla <- unique_phyla[!is.na(unique_phyla)]

carbom <- subset_taxa(carbom, (phylum %in% unique_phyla))


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
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.2) > 0, TRUE)

# HEATMAP with only most abundant OTUs
# plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

sampleOrder = sort(sample_names(carbom_abund))
# method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 

# HEATMAP time - genus, family, order, etc ...
pdf("gt_phylo_heatmap.pdf")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Species Diversity") 
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "genus", taxa.order = "genus", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Genus Diversity")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "family", taxa.order = "family", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Family Diversity")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "order", taxa.order = "order", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Order Diversity")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "class", taxa.order = "class", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Class Diversity")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE,
             taxa.label = "phylum", taxa.order = "phylum", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Phylum Diversity")
dev.off()


# HEATMAP time - cohorts
pdf("gt_phylo_heatmap_cohorts.pdf")
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  facet_grid(~ cohort, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
dev.off()

######################


# # BAR PLOT
# 
# #pdf("gt_phylo_barplot.pdf")
# # BAR GRAPH - all samples
# plot_bar(carbom, fill = "phylum") +
#   geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
#   theme(axis.text.x = element_blank())
# dev.off()
# 
# # BAR GRAPH - by time point
# plot_bar(carbom, fill = "phylum") +
#   geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
#   facet_grid(~date,scales="free_x") +
#   theme(axis.text.x = element_blank())
# pdf("gt_phylo_barplot.pdf")
# # BAR GRAPH - by time point - class
# plot_bar(carbom.abund, fill = "phylum") +
#   geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
#   #facet_grid(~date,scales="free_x") +
#   theme(axis.text.x = element_blank())
# dev.off()


######################

# DIVERSITY 

pdf("gt_phylo_diversity.pdf")
plot_richness(carbom, measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), x="cohort", color="date")
dev.off()


######################

# ORDINATION 

carbom.ord <- ordinate(carbom, "NMDS", "bray")

pdf("gt_phylo_ordination.pdf")
plot_ordination(carbom, carbom.ord, type="samples", color="date", #shape= "cohort", 
                title="gOTUs") + 
  geom_point(size=2) +
  facet_wrap(~cohort)
dev.off()

######################

# NETWORK ANALYSIS 

plot_net(carbom, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="phylum", point_label="class",
         title = "taxa network - Bray-Curtis distance")

# This is quite confusing. Let us make it more simple by using only major OTUs
plot_net(carbom_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="phylum", point_label="class",
         title = "taxa network - Bray-Curtis distance") 

pdf("gt_phylo_network.pdf")
ig = make_network(carbom_abund, type = "samples", distance = "bray", max.dist = 0.3)
plot_network(ig, carbom_abund, color = "date", shape = "cohort", line_weight = 0.3, 
             label = NULL, title = "sample network - Bray-Curtis distance")
dev.off()




######################################################################
######################################################################


# CO-HOUSING communty structure test: 

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

# we need to keep only record of pigs that were neer relocated. 
mdat2 <- setDT(mdat1)[,if(.N ==1) .SD,by=pig]


######################################################################

# SAMPLES 

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

######################################################################

# create phyloseq object

OTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)

###################################


# Control 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2","t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
ctrl <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date)
######################
# DScour 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2","t8")))
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
dscour <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date)
######################
# ColiGuard 
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2","t8")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
colig <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date)
######################
# Neomycin
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2","t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "NMDS", "bray")
neo <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date)
######################

pdf("gt_cohousing.pdf")
ggarrange(
  ctrl, neo, nrow=2, labels=c("A","B")
)
ggarrange(
  dscour, colig, nrow=2, labels=c("A","B")
)
dev.off()

