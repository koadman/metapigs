
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
library(dunn.test)
library(FSA)

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
gtdbtk_bins <- read_csv("gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

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

# load metadata (for pen info)

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


# merge pen info : 

NROW(no_reps_all)
no_reps_all <- dplyr::left_join(no_reps_all,mdat2)
NROW(no_reps_all)
head(no_reps_all)


######################################################################

# merge weight info : 

no_reps_all <- left_join(no_reps_all, weights)
NROW(no_reps_all)

######################################################################



# merge breed and bday info : 

# join breed and bday info 
NROW(no_reps_all)
no_reps_all <- left_join(no_reps_all, breed_bday)
NROW(no_reps_all)

######################################################################

# merge bins info to gtdbtk assignment info :  

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
  dplyr::select(sample,pig,date,cohort,breed,birth_day,nurse_mother,mother,weight_category,pen) %>%
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

# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 100 samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.03) > 40, TRUE)



# ORDINATION 
carbom_gt <- carbom
carbom.ord_gt <- ordinate(carbom_gt, "NMDS", "bray")

gt_ordination_plot <- plot_ordination(carbom_gt, carbom.ord_gt, type="samples", color="date") + 
  geom_point(size=1) +
  #facet_wrap(~cohort) +
  theme_bw()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 1))

pdf("gt_phylo_ordination.pdf")
gt_ordination_plot+
  facet_wrap(~cohort)+
  theme(legend.position="top")
dev.off()


########################################################################################

# NETWORK ANALYSIS 

carbom_abund_gt <- carbom_abund

ig = make_network(carbom_abund_gt, type = "samples", distance = "bray", max.dist = 0.3)
gt_network_plot <- plot_network(ig, carbom_abund_gt, color = "date", shape = "cohort", line_weight = 0.3, 
                                label = NULL, point_size = 1)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(nrow = 1))+
  guides(size = "legend", colour = "none")

pdf("gt_phylo_network.pdf")
gt_network_plot
dev.off()


########################################################################################


# HEATMAP


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

#normalization and sum of same species all in one 
z <- df %>%
  filter(!cohort=="Mothers") %>%
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

pdf("gt_phylo_diversity.pdf")
plot_richness(carbom, measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), x="cohort", color="date")
dev.off()


########################################################################################


# ORDINATION  - effect of weight 


# t0
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t2
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t4
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t6
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t6")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t6 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t8
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t8 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t10
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t10")))
carbom <- subset_samples(carbom, (weight_category %in% c("under","over")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t10 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight_category") +
  geom_point(size=2) +
  theme_bw() +
  ggtitle(paste0(sample_data(carbom)$date)[1])

weight_effect_plot <- ggarrange(p_t0,p_t2,p_t4,p_t6,p_t8,p_t10, common.legend = TRUE)

pdf("gt_phylo_ordination_weight.pdf")
weight_effect_plot
dev.off()




# # function to test correlation with Dunn.test
# myf_weight <- function(df) {
#   df1 <- df
#   weight_category <- df$weight_category
#   birth_day <- df$birth_day
#   x <- dunnTest(Axis.1~as.factor(weight_category),data=df1,method = "bonferroni")
#   x <- as.data.frame(x$res)
#   x$axis = "Axis.1"
#   y <- dunnTest(Axis.2~as.factor(weight_category),data=df1,method = "bonferroni")
#   y <- as.data.frame(y$res)
#   y$axis = "Axis.2"
#   xy <- rbind(x,y)
#   xy$comparison <- paste0(weight_category[1],"_",birth_day[1])
#   return(xy)
# }



weight_stats <- rbind(myf_weight(p_t0$data),
                      myf_weight(p_t2$data),
                      myf_weight(p_t4$data),
                      myf_weight(p_t6$data),
                      myf_weight(p_t8$data),
                      myf_weight(p_t10$data))



sink(file = "gt_phylo_ordination_weight.txt", 
     append = FALSE, type = c("output"))
paste0("Effect of weight on PCA - origin of data: phyloseq ordination data ")
paste0("PCoA  ------ Dunn (1964) Kruskal-Wallis multiple comparison
       p-values adjusted with the Bonferroni method.")
weight_stats
sink()


########################################################################################

# ORDINATION - effect of age 


######################

# age effect - all


# t0
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t2
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t4
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])


###########

dummy_all <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed")
leg_all <- get_legend(dummy_all)
all_age_and_breed <- ggarrange(p_t0,p_t2,p_t4,leg_all)


# plotted below togetehr with subset 


sink(file = "gt_phylo_ordination_age.txt",
     append = FALSE, type = c("output"))
paste0("Effect of age on PCA - ALL birth_day categories - origin of data: phyloseq ordination data ")
paste0("PC1 and PC2 respectively ------ Dunn (1964) Kruskal-Wallis multiple comparison
       p-values adjusted with the Bonferroni method.")
paste0("t0")
dunnTest(Axis.1~birth_day,data=as.data.frame(p_t0$data),method = "bonferroni")
dunnTest(Axis.2~birth_day,data=as.data.frame(p_t0$data),method = "bonferroni")
paste0("t2")
dunnTest(Axis.1~birth_day,data=as.data.frame(p_t2$data),method = "bonferroni")
dunnTest(Axis.2~birth_day,data=as.data.frame(p_t2$data),method = "bonferroni")
paste0("t4")
dunnTest(Axis.1~birth_day,data=as.data.frame(p_t4$data),method = "bonferroni")
dunnTest(Axis.2~birth_day,data=as.data.frame(p_t4$data),method = "bonferroni")
sink()

######################


# age effect - subset


# t0
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-06","2017-01-08")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t2
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-06","2017-01-08")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t4
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
carbom <- subset_samples(carbom, (birth_day %in% c("2017-01-06","2017-01-08")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])


###########

dummy_subset <- plot_ordination(carbom, carbom.ord, type="samples", color="birth_day", shape ="breed")
leg_subset <- get_legend(dummy_subset)
all_subset <- ggarrange(p_t0,p_t2,p_t4,leg_subset)


sink(file = "gt_phylo_ordination_age.txt", 
     append = TRUE, type = c("output"))
paste0("Effect of age on PCA - two birth_day categories - origin of data: phyloseq ordination data ")
paste0("PC1 and PC2 respectively ------ Dunn (1964) Kruskal-Wallis multiple comparison
       p-values adjusted with the Bonferroni method.")
paste0("t0")
dunnTest(Axis.1~birth_day,data=as.data.frame(p_t0$data),method = "bonferroni")
dunnTest(Axis.2~birth_day,data=as.data.frame(p_t0$data),method = "bonferroni")
paste0("t2")
dunnTest(Axis.1~birth_day,data=as.data.frame(p_t2$data),method = "bonferroni")
dunnTest(Axis.2~birth_day,data=as.data.frame(p_t2$data),method = "bonferroni")
paste0("t4")
dunnTest(Axis.1~birth_day,data=as.data.frame(p_t4$data),method = "bonferroni")
dunnTest(Axis.2~birth_day,data=as.data.frame(p_t4$data),method = "bonferroni")
sink()

###########

pdf("gt_phylo_ordination_age.pdf")
all_age_and_breed
all_subset
dev.off()



########################################################################################



# BREED effect - all


# t0
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t2
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])
# t4
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
p_t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="breed", shape ="breed") + 
  geom_point(size=2)  +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0(sample_data(carbom)$date)[1])


###########

dummy_all <- plot_ordination(carbom, carbom.ord, type="samples", color="breed", shape ="breed")
leg_all <- get_legend(dummy_all)
all_breed <- ggarrange(p_t0,p_t2,p_t4,leg_all)


sink(file = "gt_phylo_ordination_breed.txt", 
     append = FALSE, type = c("output"))
paste0("Effect of breed on PCA - ALL breed categories - origin of data: phyloseq ordination data ")
paste0("PC1 and PC2 respectively ------ Dunn (1964) Kruskal-Wallis multiple comparison
       p-values adjusted with the Bonferroni method.")
paste0("t0")
dunnTest(Axis.1~breed,data=as.data.frame(p_t0$data),method = "bonferroni")
dunnTest(Axis.2~breed,data=as.data.frame(p_t0$data),method = "bonferroni")
paste0("t2")
dunnTest(Axis.1~breed,data=as.data.frame(p_t2$data),method = "bonferroni")
dunnTest(Axis.2~breed,data=as.data.frame(p_t2$data),method = "bonferroni")
paste0("t4")
dunnTest(Axis.1~breed,data=as.data.frame(p_t4$data),method = "bonferroni")
dunnTest(Axis.2~breed,data=as.data.frame(p_t4$data),method = "bonferroni")
sink()



pdf("gt_phylo_ordination_breed.pdf")
all_breed
dev.off()



########################################################################################


# CO-HOUSING community structure test: 



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
carbom.ord <- ordinate(carbom, "PCoA", "bray")
ctrlt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
ctrlt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
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
carbom.ord <- ordinate(carbom, "PCoA", "bray")
dscourt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("DScour")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
dscourt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
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
carbom.ord <- ordinate(carbom, "PCoA", "bray")
coligt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
coligt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
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
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neot0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neot8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
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
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoDt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoDt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
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
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoCt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
######################
carbom <- phyloseq(OTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
# NORMALIZATION 
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoCt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  ggtitle(paste0(sample_data(carbom)$date[1],"_",sample_data(carbom)$cohort[1]))
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

sink(file = "gt_phylo_ordination_cohousing.txt", 
     append = FALSE, type = c("output"))
paste0("Effect of cohousing on PCA - origin of data: phyloseq ordination data ")
paste0("PC1 and PC2 respectively ------ Dunn (1964) Kruskal-Wallis multiple comparison
       p-values adjusted with the Bonferroni method.")
paste0("Control")
dunnTest(Axis.1~pen,data=as.data.frame(ctrlt0$data),method = "bonferroni")
dunnTest(Axis.2~pen,data=as.data.frame(ctrlt8$data),method = "bonferroni")
paste0("DScour")
dunnTest(Axis.1~pen,data=as.data.frame(dscourt0$data),method = "bonferroni")
dunnTest(Axis.2~pen,data=as.data.frame(dscourt8$data),method = "bonferroni")
paste0("ColiGuard")
dunnTest(Axis.1~pen,data=as.data.frame(coligt0$data),method = "bonferroni")
dunnTest(Axis.2~pen,data=as.data.frame(coligt8$data),method = "bonferroni")
paste0("Neomycin")
dunnTest(Axis.1~pen,data=as.data.frame(neot0$data),method = "bonferroni")
dunnTest(Axis.2~pen,data=as.data.frame(neot8$data),method = "bonferroni")
paste0("NeoD")
dunnTest(Axis.1~pen,data=as.data.frame(neoDt0$data),method = "bonferroni")
dunnTest(Axis.2~pen,data=as.data.frame(neodt8$data),method = "bonferroni")
paste0("NeoC")
dunnTest(Axis.1~pen,data=as.data.frame(neoCt0$data),method = "bonferroni")
dunnTest(Axis.2~pen,data=as.data.frame(neoCt8$data),method = "bonferroni")
sink()


########################################################################################
########################################################################################
########################################################################################
########################################################################################

############################################################################################################
############################################################################################################
############################################################################################################

carbom_ProbioticCheck <- phyloseq(OTU,TAX,samples)

############################################################################################################

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


############################################################################################################

# NORMALIZATION 

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom_ProbioticCheck))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_ProbioticCheck = transform_sample_counts(carbom_ProbioticCheck, standf)
sample_variables(carbom_ProbioticCheck)

############################################################################################################

# PLOT

######################

# HEATMAP

sampleOrder = sort(sample_names(carbom_ProbioticCheck))

# HEATMAP time - genus, family, order, etc ...
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


# probiotic strains check without phyloseq (so  lib size nornalization done before plotting)


# lib size normalisation

#################################
# STEP 1.

# normalization for library size 
df1 <- df %>%
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
getwd()


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
############################################################################################################

