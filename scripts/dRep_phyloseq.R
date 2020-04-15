
library(readr)
library(splitstackshape)
library(tidyr)
library(phyloseq)
library(dplyr)
library(ggplot2)

library(phyloseq)
library(ggpubr)
library(pheatmap)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(cluster)
library(circlize)




######################################################################

# merge info 


df <- no_reps_all %>% left_join(C1, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df)[colnames(df) == 'secondary_cluster.x'] <- 'gOTU'

NROW(unique(df$gOTU))
NROW(df)
df$secondary_cluster.y <- NULL
head(df)

######################################################################

# TAXA


taxa_mat <- df %>%
  group_by(gOTU) %>%
  slice(1) %>%
  dplyr::select(gOTU)

taxa_mat <- cSplit(taxa_mat, "gOTU", sep="_")
taxa_mat$gOTU <- paste0(taxa_mat$gOTU_1,"_",taxa_mat$gOTU_2)
taxa_mat$gOTU_2 <- taxa_mat$gOTU
NROW(taxa_mat)
NROW(unique(taxa_mat$gOTU))

taxa_mat_df <- as.data.frame(taxa_mat)
# to matrix 
taxa_mat <- taxa_mat_df
rownames(taxa_mat) <- taxa_mat[,3]
taxa_mat[,3] <- NULL
taxa_mat <- as.matrix(taxa_mat)
# ready 

NROW(unique(rownames(taxa_mat)))
head(taxa_mat_df)

######################################################################

# gOTU  
df
# columns to be kept 
keep <- c("cohort","pig","bin","date","value","gOTU")
df1 <- df[ , (names(df) %in% keep)]

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
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.2) > 0, TRUE)

# HEATMAP with only most abundant OTUs
# plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

sampleOrder = sort(sample_names(carbom_abund))
# method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 

# HEATMAP time - genus, family, order, etc ...
pdf("dRep_phylo_heatmap.pdf")
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "gOTU_1", taxa.order = "gOTU_1", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "primary clusters") 
plot_heatmap(carbom_abund, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "gOTU_2", taxa.order = "gOTU_2", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "secondary clusters") 
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

pdf("dRep_phylo_diversity.pdf")
plot_richness(carbom, measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), x="cohort", color="date")
dev.off()


######################

# ORDINATION 

carbom.ord <- ordinate(carbom, "NMDS", "bray")

pdf("dRep_phylo_ordination.pdf")
plot_ordination(carbom, carbom.ord, type="samples", color="date", #shape= "cohort", 
                title="gOTUs") + 
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

pdf("dRep_cohousing.pdf")
ggarrange(
  ctrl, neo, nrow=2, labels=c("A","B")
)
ggarrange(
  dscour, colig, nrow=2, labels=c("A","B")
)
dev.off()

########################################################################################




test <- cSplit(df2, "sample", sep="_")
colnames(test) <- c("gOTU","value","date","pig")

NROW(test)
test1 <- test

x <- test1 %>%
  filter(pig=="14159") %>%
  filter(date=="t2")
sum(x$value)

test2 <- test1 %>%
  group_by(pig,date) %>%
  dplyr::mutate(value=value/sum(value))

x <- test2 %>%
  filter(pig=="14159") %>%
  filter(date=="t2")
sum(x$value)

# to sum up 
test3 <- test2 %>%
  group_by(pig,date,gOTU) %>%
  dplyr::summarise(value=sum(value))

x <- test3 %>%
  filter(pig=="14159") %>%
  filter(date=="t2")
sum(x$value)


# add pen info
NROW(test3)
test4 <- dplyr::left_join(test3,mdat2)
NROW(test4)

test4$sample <- paste0(test4$pen,"_",test4$pig)
head(test4)

NROW(unique(test4$gOTU))

# you have to remove (now that you normalized it's fine) all the samples without 
# a pen description
NROW(which(is.na(test4)))

NROW(test4)
test5 <- na.omit(test4) %>%
  filter(date=="t1"|date=="t3"|date=="t5"|date=="t7")
NROW(test5)

unique(test5$gOTU)

myg <- test5 %>%
  group_by(gOTU) %>%
  filter(n()>80)  %>%
  filter(n()<100)
myg <- unique(myg$gOTU)
NROW(myg)

pdf("plotssss.pdf")
for (i in myg) {
  
  toplot <- test5 %>%
    filter(gOTU==i) %>%
    dplyr::select(date,sample,value)  %>%
    pivot_wider(names_from = date, values_from = value, values_fill = list(value = 0)) 
  
  if (NROW(toplot) > 6) {
    
    toplot[1] <- NULL
    toplot <- as.data.frame(toplot)
    
    
    m <- as.matrix(toplot[, -1])
    rownames(m) <- toplot$sample
    colnames(m) <- sort(colnames(m))
    rownames(m) <- sort(rownames(m))
   
    # to binary
    # m <- ifelse(m>0, 1, 0)
    Heatmap(m, name = "mat", row_split = factor(substr(rownames(m), start = 1, stop = 2)),
            #row_km = 3, 
            row_gap = unit(2, "mm"))
    #heatmap(m, Rowv = NA, Colv = NA, main = paste0("secondary_cluster: ",i), na.rm = FALSE)
  }
  
}
dev.off()


Heatmap(m, name = "mat", row_split = factor(substr(rownames(m), start = 1, stop = 2)),
        #row_km = 3, 
        row_gap = unit(2, "mm"))
