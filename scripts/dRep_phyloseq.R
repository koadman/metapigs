
library(readr)
library(splitstackshape)
library(tidyr)
library(phyloseq)
library(dplyr)
library(ggplot2)

library(phyloseq)
library(ggpubr)
library(pheatmap)

#library(ComplexHeatmap)
library(cluster)
library(circlize)
library(readxl)
library(data.table)
library(FSA)
library(openxlsx)



setwd("~/Desktop/metapigs_dry/dRep/")
basedir = "~/Desktop/metapigs_dry/"


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

# ###### run this part if you want to only retain piglets that 
# # where present throughput trial (not euthanised) 
# 
# to_keep <- no_reps_all %>%
#   filter(date=="t8") %>%
#   dplyr::select(pig) %>%
#   distinct()
# 
# to_keep <- as.character(to_keep$pig)
# 
# no_reps_all <- subset(no_reps_all, (pig %in% to_keep))
# NROW(unique(no_reps_all$pig))
# 
# ############################################################

######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(basedir,"gtdbtk/gtdbtk_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################
######################################################################


# upload metadata for pen info

mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)

mdat2 <- mdat %>%
  dplyr::filter(!Cohort=="Sows")  %>%
  dplyr::filter(!`*collection_date`=="2017-01-31"|
                  `*collection_date`=="2017-02-01"|
                  `*collection_date`=="2017-02-03") %>%
  dplyr::select(isolation_source,PigPen)

colnames(mdat2) <- c("pig","pen")
mdat2 <- as.data.frame(mdat2)

mdat2 <- mdat2 %>%
  group_by(pig) %>%
  distinct()
NROW(mdat2)

mdat2$pen <- gsub("nan",NA,mdat2$pen)
mdat2 <- na.omit(mdat2)

# we need to keep only record of pigs that were not relocated. 
mdat2 <- setDT(mdat2)[,if(.N ==1) .SD,by=pig]


######################################################################


# upload breed and bday info 

suppl_piglets_details_mothers_weight <- read_excel("~/Desktop/metapigs_dry/suppl_piglets_details_mothers&weight.xlsx")

# select cols of interest
breed_bday <- suppl_piglets_details_mothers_weight %>%
  dplyr::select(TATTOO,BIRTH_DAY,...8,`Nursing Dam`,STIGDAM)

# rename columns
colnames(breed_bday) <- c("pig","birth_day","breed","nurse_mother","mother")

breed_bday$birth_day <- as.character(breed_bday$birth_day)

# clean names
breed_bday$pig <- gsub("G","", breed_bday$pig)
breed_bday$pig <- gsub("T","", breed_bday$pig)

breed_bday <- as.data.frame(breed_bday)

######################################################################

# upload weight info 


weights <- read_csv(paste0(basedir,"weights.csv"), 
                    col_types = cols(Pig = col_character(), 
                                     Room = col_character()))
colnames(weights) <- c("room","pen","pig","t0","t2","t4","t6","t8")


weights_final <- read_csv(paste0(basedir,"weights_final.csv"), 
                          col_types = cols(Pig = col_character(), 
                                           Room = col_character()))
colnames(weights_final) <- c("room","pen","pig","date","weight")
weights_final$date <- gsub("6-Mar","t10",weights_final$date)
weights_final$date <- gsub("7-Mar","t10",weights_final$date)
weights_final$date <- gsub("8-Mar","t10",weights_final$date)
weights_final$date <- gsub("9-Mar","t10",weights_final$date)
weights_final <- weights_final %>%
  dplyr::select(pig,date,weight) %>%
  filter(!date=="10-Mar") # as it's NA

weights <- weights %>%
  dplyr::select(pig,t0,t2,t4,t6,t8) %>%
  pivot_longer(., cols = c(t0,t2,t4,t6,t8), names_to = "date", values_to = "weight")
weights <- as.data.frame(weights)

weights <- rbind(weights,weights_final)
NROW(weights)

# merge bday info : 
bday <- breed_bday %>%
  dplyr::select(pig,birth_day)
weights <- left_join(weights, bday)
NROW(weights)
unique(weights$birth_day)

weights_rest <- weights %>% 
  filter(!birth_day=="2017-01-06") %>%
  group_by(birth_day,date) %>%
  mutate(weight_category=cut(weight, breaks=c(summary(weight)[1], summary(weight)[2], summary(weight)[5], summary(weight)[6]), 
                             labels=c("under","normal","over"))) 

weights_rest<- as.data.frame(weights_rest)

# quickly visualize weight category distribution
ggplot(weights_rest,aes(x=date,fill=weight_category)) +
  geom_bar() +
  facet_wrap(~birth_day)

weights6 <- weights %>% 
  filter(birth_day=="2017-01-06")
weights6$weight_category <- NA

weights <- rbind(weights_rest,weights6) %>%
  dplyr::select(pig,date,weight_category)
head(weights)


######################################################################

# merge bins info to gtdbtk assignment info :  

NROW(gtdbtk_bins)
NROW(no_reps_all)
head(gtdbtk_bins)
head(no_reps_all)
df0 <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df0)[colnames(df0) == 'node'] <- 'gOTU'
df0$gOTU <- paste0("gOTU_",df0$gOTU)

NROW(unique(df0$gOTU))
NROW(df0)

######################################################################

# merge all otehr info: 
# add pen info (mdat2), breed and bday info (breed_bday) and weight info (weights)

# add breed and bday info (breed_bday)
df0 <- left_join(df0,breed_bday)
NROW(df0)

# add pen info (mdat2)
df0 <- left_join(df0,mdat2)
NROW(df0)

# add weight info (weights)
df0 <- left_join(df0,weights)
NROW(df0)

###########################################################################################

# create workbook to add stats 

wb <- createWorkbook()

###########################################################################################

# TAXA


taxa_mat <- df0 %>%
  dplyr::select(secondary_cluster,species,genus,family,order,class,phylum,domain) %>%
  group_by(secondary_cluster) %>%
  slice(1)

NROW(taxa_mat)
NROW(unique(taxa_mat$secondary_cluster))

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

# species  

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","secondary_cluster")
df1 <- df0[ , (names(df0) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"

NROW(df1)
# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,date,secondary_cluster) %>%
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

# check whether this list is empty(no NAs)
check_DF <- df3[rowSums(is.na(df3)) > 0,]
NROW(check_DF)


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
  dplyr::select(sample,pig,date,cohort,pen,birth_day,breed,weight_category) %>%
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

carbom_dRep <- phyloseq(OTU,TAX,samples)

############################################################################################################

# SUBSETTING phyloseq obejct

carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9")))


###########################################################################################

# NORMALIZATION 

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)

############################################################################################################

# PLOT

######################


# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 40 samples 
carbom_dRep_abund <- filter_taxa(carbom_dRep, function(x) sum(x > total*0.03) > 40, TRUE)


# ORDINATION 


carbom_dRep_abund.ord <- ordinate(carbom_dRep_abund, "NMDS", "bray")

dRep_ordination_plot <- plot_ordination(carbom_dRep_abund, carbom_dRep_abund.ord, type="samples", color="date") + 
  geom_point(size=1) +
  #facet_wrap(~cohort) +
  theme_bw() +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 1))

pdf("dRep_phylo_ordination.pdf")
dRep_ordination_plot+
  facet_wrap(~cohort)+
  theme(legend.position="top")
dev.off()

######################

# NETWORK ANALYSIS 
  
ig = make_network(carbom_dRep_abund, type = "samples", distance = "bray", max.dist = 0.3)
dRep_network_plot <- plot_network(ig, carbom_dRep_abund, color = "date", shape = "cohort", line_weight = 0.3, 
                                label = NULL, point_size = 1)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(nrow = 1))+
  guides(size = "legend", colour = "none")

pdf("dRep_phylo_network.pdf")
dRep_network_plot
dev.off()


######################

# HEATMAP - with only most abundant OTUs

sampleOrder = sort(sample_names(carbom_dRep_abund))
# method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 

# HEATMAP time - genus, family, order, etc ...
pdf("dRep_phylo_heatmap.pdf")
plot_heatmap(carbom_dRep_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "gOTU_2", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster") 
plot_heatmap(carbom_dRep_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Species)")
plot_heatmap(carbom_dRep_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "genus", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Genus)")
plot_heatmap(carbom_dRep_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "family", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Family)")
plot_heatmap(carbom_dRep_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "order", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Order)")
plot_heatmap(carbom_dRep_abund, method = "MDS", distance="unifrac",weighted=TRUE,
             taxa.label = "class", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Class)")
plot_heatmap(carbom_dRep_abund, method = "MDS", distance="unifrac",weighted=TRUE,
             taxa.label = "phylum", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Microbe Diversity by secondary cluster (Phylum)")
dev.off()


######################

# DIVERSITY 

# here rarefaction is applied

carbom <- phyloseq(OTU,TAX,samples)

carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9")))

carbom_rarefied = rarefy_even_depth(carbom, replace=TRUE, rngseed = 42)

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom_rarefied))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_rarefied = transform_sample_counts(carbom_rarefied, standf)
sample_variables(carbom_rarefied)


dRep_diversity <- plot_richness(carbom_rarefied, measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), 
                              x="cohort", color="date") +
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

dRep_diversity_small <- plot_richness(carbom_rarefied, measures=c("Chao1","Shannon","InvSimpson"), 
                                    color="date") +
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top",
        axis.text.x = element_blank())

pdf("dRep_phylo_diversity.pdf")
dRep_diversity
dev.off()

######################################################################
######################################################################


# ordination - CO-HOUSING : 


# Control 
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t0")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("Control")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
ctrlt0 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t8")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("Control")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
ctrlt8 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# DScour 
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t0")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("DScour")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
dscourt0 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t8")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("DScour")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
dscourt8 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# ColiGuard 
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t0")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("ColiGuard")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
coligt0 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t8")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("ColiGuard")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
coligt8 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# Neomycin
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t0")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("Neomycin")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
neot0 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t8")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("Neomycin")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
neot8 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoD 
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t0")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("NeoD")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
neoDt0 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t8")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("NeoD")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
neoDt8 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoC 
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t0")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("NeoC")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
neoCt0 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom_dRep <- phyloseq(OTU,TAX,samples)
carbom_dRep <- subset_samples(carbom_dRep, (date %in% c("t8")))
carbom_dRep <- subset_samples(carbom_dRep, (cohort %in% c("NeoC")))
# NORMALIZATION 
total = median(sample_sums(carbom_dRep))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_dRep = transform_sample_counts(carbom_dRep, standf)
sample_variables(carbom_dRep)
carbom_dRep.ord <- ordinate(carbom_dRep, "PCoA", "bray")
neoCt8 <- plot_ordination(carbom_dRep, carbom_dRep.ord, type="samples", color="pen", shape="pen",
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

pdf("dRep_phylo_ordination_cohousing.pdf")
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


# function to test correlation with Dunn.test

pen_function <- function(myd){
  
  x_axis1 <- dunnTest(myd$Axis.1~pen,data=myd,method = "bonferroni")
  x_axis1 <- as.data.frame(x_axis1$res)
  x_axis1$axis = "Axis.1"
  
  x_axis2 <- dunnTest(myd$Axis.2~pen,data=myd,method = "bonferroni")
  x_axis2 <- as.data.frame(x_axis2$res)
  x_axis2$axis = "Axis.2"
  
  x <- rbind(x_axis1,x_axis2)
  
}

out1 <- pen_function(ctrlt0$data)
out1$type = "ctrlt0"
out2 <- pen_function(ctrlt8$data)
out2$type = "ctrlt8"
out3 <- pen_function(dscourt0$data)
out3$type = "dscourt0"
out4 <- pen_function(dscourt8$data)
out4$type = "dscourt8"
out5 <- pen_function(coligt0$data)
out5$type = "coligt0"
out6 <- pen_function(coligt8$data)
out6$type = "coligt8"
out7 <- pen_function(neot0$data)
out7$type = "neot0"
out8 <- pen_function(neot8$data)
out8$type = "neot8"
out9 <- pen_function(neoDt0$data)
out9$type = "neoDt0"
out10 <- pen_function(neoDt8$data)
out10$type = "neoDt8"
out11 <- pen_function(neoCt0$data)
out11$type = "neoCt0"
out12 <- pen_function(neoCt8$data)
out12$type = "neoCt8"


pen_stats <- rbind(out1,out2,out3,out4,out5,out6,
      out7,out8,out9,out10,out11,out12)
pen_stats$test = "Dunn.test"
pen_stats$correction = "Bonferroni"


addWorksheet(wb, "cohousing_dRep")
writeData(wb, sheet = "cohousing_dRep", pen_stats, rowNames = FALSE)

########################################################################################


# save stats in workbook
saveWorkbook(wb, paste0(basedir,"stats.xlsx"), overwrite=TRUE)


