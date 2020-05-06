
library(readr)
library(splitstackshape)
library(tidyr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(EnvStats)
library(ggpubr)
library(readxl)
library(data.table)
library(ape)
library(scales)
library(FSA)
library(openxlsx)

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
# gt_phylo_heatmap_ProbioticCheck.pdf
# gt_ProbioticCheck.pdf
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

# merge all other info: 
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

# gOTU_mat

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","gOTU")
df1 <- df0[ , (names(df0) %in% keep)]

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

# check whether this list is empty(no NAs)
check_DF <- df3[rowSums(is.na(df3)) > 0,]
NROW(check_DF)


# to matrix 
gOTU_mat <- as.data.frame(df3)
rownames(gOTU_mat) <- gOTU_mat[,1]
gOTU_mat[,1] <- NULL
gOTU_mat <- as.matrix(gOTU_mat)
mode(gOTU_mat) <- "integer"
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

unique(sample_df$breed)

sample_df$breed <- gsub("Landrace x Cross bred [(]LW x D[])]","LxLWD", sample_df$breed)
sample_df$breed <- gsub("Duroc x Landrace","DxL", sample_df$breed)
sample_df$breed <- gsub("Duroc x Large white","DxLW", sample_df$breed)
sample_df$breed <- gsub("Large white x Duroc","LWxD", sample_df$breed)

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

gOTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)


############################################################################################################

# PLOT

######################


# ORDINATION 

# NORMALIZATION BY MEDIAN SEQUENCING DEPTH
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9")))
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 40 samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.03) > 40, TRUE)

carbom_abund.ord <- ordinate(carbom_abund, "NMDS", "bray")

gt_ordination_plot <- plot_ordination(carbom_abund, carbom_abund.ord, type="samples", color="date") + 
  geom_point(size=1) +
  theme_bw()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

pdf("gt_phylo_ordination.pdf")
gt_ordination_plot +
  facet_wrap(~cohort)
dev.off()


########################################################################################

# NETWORK ANALYSIS 

# NORMALIZATION BY MEDIAN SEQUENCING DEPTH
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9")))
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 40 samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.03) > 40, TRUE)

ig = make_network(carbom_abund, type = "samples", distance = "bray", max.dist = 0.35)
gt_network_plot <- plot_network(ig, carbom_abund, color = "date", shape = "cohort", line_weight = 0.3, 
                                label = NULL, point_size = 1)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(nrow = 1))+
  guides(size = "legend", colour = "none")

pdf("gt_phylo_network.pdf")
gt_network_plot
dev.off()


########################################################################################


# HEATMAP

# NORMALIZATION BY RAREFACTION
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# keep only very abundant OTUs: more than 5 counts per sample, in at least 1/4 samples 
carbom_abund <- filter_taxa(carbom, 
                            function(x) 
                              sum(x > 5) > (NROW(sample_data(carbom))/4), 
                            TRUE)

random_tree = rtree(ntaxa(carbom_abund), rooted=TRUE, tip.label=taxa_names(carbom_abund))
physeq1 = merge_phyloseq(carbom_abund,random_tree)


# HEATMAP time - genus, family, order, etc ...
pdf("gt_phylo_heatmap.pdf")
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Species Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "genus", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Genus Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "family", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Family Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "order", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Order Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "class", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Class Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "phylum", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Phylum Diversity") 
dev.off()


# HEATMAP time - cohorts
pdf("gt_phylo_heatmap_cohorts.pdf")
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  facet_grid(~ cohort, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
dev.off()


##############################################

# this is a small break in between the phyloseq analysis: 

# DISPLAYING MOST ABUNDANT SPECIES TIME-TREND

# this is a zoom in on the species shown in phylo_heatmap because most abundant 
# (taking gOTUs that represent at least 3% of the sample and present in at least 40 samples)

# 1. raw data is normalized by lib size and 
# 2. the mean is taken from bins assigned the same species

# 3. at this point those species are now selected and plotted using the log10

physeq1
keep_these_taxa <- as.data.frame(tax_table(physeq1))$species
keep_these_taxa <- as.character(keep_these_taxa)

# normalization and sum of same species all in one 
z <- df0 %>%
  dplyr::filter(!cohort=="Mothers") %>%
  group_by(pig,date) %>%
  mutate(norm_value=value/sum(value)) %>% # lib size normalization
  group_by(cohort,pig,date,species) %>%
  dplyr::summarize(z = mean(norm_value)) # mean of bins falling within same species

z <- subset(z, (species %in% keep_these_taxa))
NROW(unique(z$species))

z <- as.data.frame(z)

# reorder dates 
z$date  = factor(z$date, levels=c("t0",
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
z$cohort  = factor(z$cohort, levels=c("Control",
                                      "DScour", 
                                      "ColiGuard",
                                      "Neomycin",
                                      "NeoD",
                                      "NeoC"))


# split df by species
multiple_DFs <- split( z , f = z$species )

pdf("gt_zoomIN_on_phylo_heatmap.pdf")
for (single_DF in multiple_DFs) {
  
  DF <- as.data.frame(single_DF)
  this_species <- DF$species
  
  p <- ggplot(DF, aes(x=date, y=log10(z), fill=cohort)) + 
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(size=0.2)+       # keep this if you want to see how many samples behind measurement
    #geom_line() + geom_point(size=0.8)+
    theme_bw()+
    theme(legend.position="right")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "collection date",
         y = "relative abundance (log10)",
         color = "Phylum")  +
    theme(legend.title = element_text()) +
    ggtitle(this_species)
  
  print(p)
  
}
dev.off()


##############################################


# BAR PLOT


# whether you run it with 
# median sequencing depth normalization or rarefaction
# output doesn t change much except the Actinobacteriota not showing when you rarefy 
# trend is the same with either normalization method 

# NORMALIZATION BY MEDIAN SEQUENCING DEPTH
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9")))
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 40 samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.03) > 40, TRUE)

# BAR GRAPH - by time point
pdf("gt_phylo_barplot_time.pdf")
plot_bar(carbom_abund, fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
plot_bar(carbom_abund, fill = "class") +
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
dev.off()


######################

# DIVERSITY 

# NORMALIZATION BY RAREFACTION
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# keep only very abundant OTUs: 
# more than 2 counts in at least a fourth of the toal number of samples 
carbom_abund <- filter_taxa(carbom, 
                            function(x) 
                              sum(x > 2) > (NROW(sample_data(carbom))/4), 
                            TRUE)

gt_diversity_samples <- plot_richness(carbom, measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), 
                                      color="date", x="date") +
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")


pdf("gt_phylo_diversity.pdf")
gt_diversity_samples
dev.off()


######### plotting above results in a different way: 
# focus on three measures;
# whisker plots instead 

my_comparisons <- list( c("t0", "t2"), c("t2", "t4"), c("t4", "t6"),
                        c("t6", "t8"), c("t4", "t8"))

Chao1 <- gt_diversity_samples$data %>%
  filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  filter(variable=="Chao1") %>%
  ggplot(., aes(x=date, y=value, color=date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Chao1") +
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif") # Add pairwise comparisons p-value

Shannon <- gt_diversity_samples$data %>%
  filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  filter(variable=="Shannon") %>%
  ggplot(., aes(x=date, y=value, color=date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Shannon") +
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif") # Add pairwise comparisons p-value

Simpson <- gt_diversity_samples$data %>%
  filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  filter(variable=="Simpson") %>%
  ggplot(., aes(x=date, y=value, color=date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Simpson") +
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif") # Add pairwise comparisons p-value

tosave <- ggarrange(Chao1, 
                    Shannon,
                    Simpson,
                    ncol = 3, 
                    nrow=1,
                    labels=c("A","B","C"), 
                    common.legend = TRUE)

ggsave(filename = "gt_phylo_diversity_boxplot.pdf", plot = tosave)


######################################################################
######################################################################


# ordination - CO-HOUSING : 


# Control 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
ctrlt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
ctrlt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# DScour 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
dscourt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
dscourt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# ColiGuard 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
coligt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
coligt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# Neomycin
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neot0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neot8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoD 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoDt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoDt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoC 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoCt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")

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

pdf("gt_phylo_ordination_cohousing.pdf")
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


addWorksheet(wb, "cohousing_gtdb")
writeData(wb, sheet = "cohousing_gtdb", pen_stats, rowNames = FALSE)



######################################################################
######################################################################


# ORDINATION  - effect of weight 


# t0
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                            function(x) 
                              sum(x > 10) > (NROW(sample_data(carbom))/2), 
                            TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t2
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t4
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t6
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t6")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t6 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t8
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t8 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t10
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t10")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t10 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])

weight_effect_plot <- ggarrange(p_t0,p_t2,p_t4,p_t6,p_t8,p_t10, common.legend = TRUE)

pdf("gt_phylo_ordination_weight.pdf")
weight_effect_plot
dev.off()



# function to test correlation with Dunn.test

weight_function <- function(myd){
  
  x_axis1 <- dunnTest(myd$Axis.1~weight_category,data=myd,method = "bonferroni")
  x_axis1 <- as.data.frame(x_axis1$res)
  x_axis1$axis = "Axis.1"
  
  x_axis2 <- dunnTest(myd$Axis.2~weight_category,data=myd,method = "bonferroni")
  x_axis2 <- as.data.frame(x_axis2$res)
  x_axis2$axis = "Axis.2"
  
  x <- rbind(x_axis1,x_axis2)
  
}

out1 <- weight_function(p_t0$data)
out1$type = "t0"
out2 <- weight_function(p_t2$data)
out2$type = "t2"
out3 <- weight_function(p_t4$data)
out3$type = "t4"
out4 <- weight_function(p_t6$data)
out4$type = "t6"
out5 <- weight_function(p_t8$data)
out5$type = "t8"
out6 <- weight_function(p_t10$data)
out6$type = "t10"


weight_stats <- rbind(out1,out2,out3,out4,out5,out6)

weight_stats$test = "Dunn.test"
weight_stats$correction = "Bonferroni"


addWorksheet(wb, "weight_gtdb")
writeData(wb, sheet = "weight_gtdb", weight_stats, rowNames = FALSE)



######################################################################
######################################################################


# ordination - AGE (birth_day) : 

# age effect - all piglets 

# t0
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
a_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~breed)
dummy_all <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() 
# t2
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
a_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~breed)
# t4
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
a_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~breed)


###########

leg_all <- get_legend(dummy_all)
all_age_and_breed <- ggarrange(a_t0,a_t2,a_t4,leg_all)


######################

# age effect - subset

# t0
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-08","2017-01-11")))
carbom <- subset_samples(carbom, (breed %in% c("DxL")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
s_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~breed)
dummy_subset <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() 
# t2
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-08","2017-01-11")))
carbom <- subset_samples(carbom, (breed %in% c("DxL")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
s_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~breed)
# t4
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-08","2017-01-11")))
carbom <- subset_samples(carbom, (breed %in% c("DxL")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
s_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~breed)


###########

leg_subset <- get_legend(dummy_subset)
all_subset <- ggarrange(s_t0,s_t2,s_t4,leg_subset)

###########

pdf("gt_phylo_ordination_age.pdf")
all_age_and_breed
all_subset
dev.off()

###########

# stats

# function to test correlation with Dunn.test

age_function <- function(myd){
    
    x_axis1 <- dunnTest(myd$Axis.1~birth_day,data=myd,method = "bonferroni")
    x_axis1 <- as.data.frame(x_axis1$res)
    x_axis1$axis = "Axis.1"
    
    x_axis2 <- dunnTest(myd$Axis.2~birth_day,data=myd,method = "bonferroni")
    x_axis2 <- as.data.frame(x_axis2$res)
    x_axis2$axis = "Axis.2"
    
    x <- rbind(x_axis1,x_axis2) 
    x$group <- as.character(a_t2$data$date[1])
    
    x <- as.data.frame(x)
}

out1 <- age_function(a_t0$data)
out1$group_analyzed = "all"
out2 <- age_function(a_t2$data)
out2$group_analyzed = "all"
out3 <- age_function(a_t4$data)
out3$group_analyzed = "all"
out4 <- age_function(s_t0$data)
out4$group_analyzed = "subset"
out5 <- age_function(s_t2$data)
out5$group_analyzed = "subset"
out6 <- age_function(s_t4$data)
out6$group_analyzed = "subset"


age_stats <- rbind(out1,out2,out3,out4,out5,out6)

age_stats$test = "Dunn.test"
age_stats$correction = "Bonferroni"


addWorksheet(wb, "age_gtdb")
writeData(wb, sheet = "age_gtdb", age_stats, rowNames = FALSE)



######################################################################
######################################################################


# ordination - BREED : 


#######

# t0
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~birth_day)
dummy_all <- plot_ordination(carbom, carbom.ord, type="samples", color="breed")
# t2
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~birth_day)
# t4
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~birth_day)

###########

leg_all <- get_legend(dummy_all)
all_breed <- ggarrange(p_t0,p_t2,p_t4,leg_all)

######################

# age effect - subset

# t0
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (breed %in% c("DxL","DxLW")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-08")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
s_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~birth_day)
dummy_subset <- plot_ordination(carbom, carbom.ord, type="samples", color="breed") + 
  geom_point(size=2)  +
  theme_bw() 
# t2
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
carbom <- subset_samples(carbom, (breed %in% c("DxL","DxLW")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-08")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
s_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~birth_day)
# t4
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
carbom <- subset_samples(carbom, (breed %in% c("DxL","DxLW")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-08")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom <- filter_taxa(carbom, 
                      function(x) 
                        sum(x > 10) > (NROW(sample_data(carbom))/2), 
                      TRUE)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
s_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])+
  facet_wrap(~birth_day)


###########

leg_subset <- get_legend(dummy_subset)
breed_subset <- ggarrange(s_t0,s_t2,s_t4,leg_subset)

###########

pdf("gt_phylo_ordination_breed.pdf")
all_breed
breed_subset
dev.off()



# function to test correlation with Dunn.test

breed_function <- function(myd){
  
  x_axis1 <- dunnTest(myd$Axis.1~breed,data=myd,method = "bonferroni")
  x_axis1 <- as.data.frame(x_axis1$res)
  x_axis1$axis = "Axis.1"
  
  x_axis2 <- dunnTest(myd$Axis.2~breed,data=myd,method = "bonferroni")
  x_axis2 <- as.data.frame(x_axis2$res)
  x_axis2$axis = "Axis.2"
  
  x <- rbind(x_axis1,x_axis2)
  
}

out1 <- breed_function(p_t0$data)
out1$group = "t0"
out1$group_analyzed = "all"
out2 <- breed_function(p_t2$data)
out2$group = "t2"
out2$group_analyzed = "all"
out3 <- breed_function(p_t4$data)
out3$group = "t4"
out3$group_analyzed = "all"
out4 <- breed_function(s_t0$data)
out4$group = "t0"
out4$group_analyzed = "subset"
out5 <- breed_function(s_t2$data)
out5$group = "t2"
out5$group_analyzed = "subset"
out6 <- breed_function(s_t4$data)
out6$group = "t4"
out6$group_analyzed = "subset"


breed_stats <- rbind(out1,out2,out3,out4,out5,out6)

breed_stats$test = "Dunn.test"
breed_stats$correction = "Bonferroni"


addWorksheet(wb, "breed_gtdb")
writeData(wb, sheet = "breed_gtdb", breed_stats, rowNames = FALSE)


age_stats_sign <- age_stats %>%
  dplyr::select(P.adj,group_analyzed,group,Comparison) %>%
  filter(P.adj<0.05) %>%
  arrange(P.adj)
  
breed_stats_sign <- breed_stats %>%
  dplyr::select(P.adj,group_analyzed,group,Comparison) %>%
  filter(P.adj<0.05) %>%
  arrange(P.adj)

sink("breed_age_significant.txt", type = "output", append = FALSE)
start_message <- " ########################## AGE - significant ########################## "
start_message
age_stats_sign
start_message <- " ########################## BREED - significant ########################## "
start_message
breed_stats_sign
sink()

############################################################################################################
############################################################################################################
############################################################################################################

# Presence of PROBIOTIC species in dataset - plotting only those 


carbom_ProbioticCheck <- phyloseq(gOTU,TAX,samples)
carbom_ProbioticCheck <- subset_samples(carbom_ProbioticCheck, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9")))
carbom_ProbioticCheck <- subset_taxa(carbom_ProbioticCheck, (species %in% c("Lactobacillus_B salivarius",
                                                                            "Lactobacillus_F plantarum",
                                                                            "Enterococcus_B faecium_B",
                                                                            "Bifidobacterium bifidum",
                                                                            "Streptococcus thermophilus",
                                                                            "Lactobacillus_F plantarum", 
                                                                            "Lactobacillus helveticus",
                                                                            "Lactobacillus delbreuckii",
                                                                            "Lactobacillus_C rhamnosus")))
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom_ProbioticCheck))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_ProbioticCheck = transform_sample_counts(carbom_ProbioticCheck, standf)
sample_variables(carbom_ProbioticCheck)

# HEATMAP

sampleOrder = sort(sample_names(carbom_ProbioticCheck))
pdf("gt_phylo_heatmap_ProbioticCheck.pdf")
# HEATMAP time - cohorts
plot_heatmap(carbom_ProbioticCheck, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  facet_grid(~ cohort, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Probiotic Species - Presence in cohorts")
dev.off()



######################################################################
######################################################################


# PROBIOTIC species check without phyloseq (so the lib size normalization is done before plotting)


# lib size normalisation

#################################
# STEP 1.

# normalization for library size 
df1 <- df0 %>%
  dplyr::group_by(pig,date) %>% # a pig sample is: pig+date
  dplyr::mutate(norm_value=value/sum(value))    
head(df1)
sort(unique(df1$species))


df2 <- df1 %>%
  filter(!date=="tM") %>%
  filter(!date=="t10") %>%
  group_by(cohort,pig,date,species) %>%
  dplyr::summarise(disp = mean(norm_value), sd = sd(norm_value), n=n()) %>%
  filter(species=="Lactobacillus_B salivarius"|
           species=="Lactobacillus_F plantarum"|
           species=="Bifidobacterium bifidum"|
           species=="Enterococcus_B faecium_B"|
           species=="Lactobacillus delbrueckii"|
           species=="Lactobacillus helveticus"|
           species=="Lactobacillus_C rhamnosus"|
           species=="Lactobacillus_F plantarum"|
           species=="Streptococcus thermophilus")

# reorder cohorts 
df2$cohort  = factor(df2$cohort, levels=c("Control", 
                                          "DScour",
                                          "ColiGuard", 
                                          "Neomycin",
                                          "NeoD",
                                          "NeoC"))

multiple_DFs <- split( df2 , f = df2$species )



pdf("gt_ProbioticCheck.pdf")
for (single_DF in multiple_DFs) {
  
  
  DF <- as.data.frame(single_DF)
  
  species <- DF$species[1]
  DF$species <- NULL
  
  print(ggplot(DF, aes(x=date,y=disp, color=cohort)) +
          #geom_bar(stat="identity",position=position_dodge())+
          geom_point(aes(x=date,y=disp, group=pig),size=2)+
          geom_line(aes(x=date,y=disp, group=pig))+
          theme_minimal()+
          stat_n_text(size = 2) +
          facet_wrap(~cohort, scales = c("fixed")) +
          ggtitle(paste0(species)))
  
  
}
dev.off()



############################################################################################################
############################################################################################################

# save stats in workbook
saveWorkbook(wb, paste0(basedir,"stats.xlsx"), overwrite=TRUE)

############################################################################################################
############################################################################################################


