
library(readr)
library(tidyverse)
library(dplyr)
library(robCompositions)
library(microbiome)
library(phyloseq)
library(DESeq2)
library(ggplot2)

######################################################################

# counts data 

no_reps_all <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/no_reps_all.csv", 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)


######################################################################

# checkM data 

# upload input 
concatenated_checkm_nearly <- read_delim("/Users/12705859/Desktop/metapigs_dry/checkm/checkm_all_nearly", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)

NROW(concatenated_checkm_nearly)
# some formatting 
concatenated_checkm_nearly$Completeness <- as.numeric(concatenated_checkm_nearly$Completeness)
concatenated_checkm_nearly$Contamination <- as.numeric(concatenated_checkm_nearly$Contamination)
# remove rows containing NAs as in the original file these rows where headers
concatenated_checkm_nearly <- na.omit(concatenated_checkm_nearly)

# filter >90 <5 bins 
concatenated_checkm_nearly <- dplyr::filter(concatenated_checkm_nearly, !grepl("Completeness",Completeness))
newdata <- subset(concatenated_checkm_nearly, Completeness >= 90 & Contamination <= 5)

# rename cols to matching colnames between dataframes to merge  
colnames(newdata)[colnames(newdata)=="pigid"] <- "pig"
colnames(newdata)[colnames(newdata)=="Bin Id"] <- "bin"
colnames(newdata)[colnames(newdata)=="Taxonomy (contained)"] <- "taxa"

newdata <- newdata %>%
  select(pig,bin,taxa)

######################################################################


# create "OTU" table: a separate "OTU" identifier for each different taxa

new <- newdata

# split taxa column into several (kingdom, phylum, etc ...) 
new$all_taxa <- new$taxa
new <- cSplit(new, "taxa", sep=";")
head(new)
colnames(new) <- c("pig","bin","all_taxa","kingdom","phylum","class","order","family","genus","species")
head(new)
NROW(new)

new <- new %>%
  group_by(all_taxa) %>%
  mutate(OTU = paste0("OTU_",group_indices())) 

taxa_mat <- new[4:11] 

NROW(taxa_mat)
NROW(unique(taxa_mat$OTU))
taxa_mat <- taxa_mat %>%
  distinct()
NROW(taxa_mat)
NROW(unique(taxa_mat$OTU))

taxa_mat_df <- as.data.frame(taxa_mat)
# to matrix 
taxa_mat <- taxa_mat_df
rownames(taxa_mat) <- taxa_mat[,8]
taxa_mat[,8] <- NULL
taxa_mat <- as.matrix(taxa_mat)
# ready 

NROW(unique(rownames(taxa_mat)))
head(taxa_mat_df)

############################################################################################################

df <- right_join(no_reps_all, new, by=c("pig","bin"))
head(df)
NROW(df)

NROW(unique(paste0(df$pig,df$date)))
NROW(unique(paste0(df$OTU)))

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","OTU")
df1 <- df[ , (names(df) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"

NROW(df1)
# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,date,OTU) %>%
  dplyr::summarise(all_bins_value = sum(value))

NROW(df2)
NROW(unique(paste0(df2$pig,df2$date)))


# assign a unique sample name 
df2$sample <- paste0(df2$pig,"_",df2$date)
# remove now pig and date (redundant)
df2 <- df2[3:5]

# long to wide 
df3 <- df2 %>%
  pivot_wider(names_from = sample, values_from = all_bins_value, values_fill = list(all_bins_value = 0)) 

# to matrix 
otu_mat <- as.data.frame(df3)
rownames(otu_mat) <- otu_mat[,1]
otu_mat[,1] <- NULL
otu_mat <- as.matrix(otu_mat)
# ready to transform to OTU table

NROW(unique(rownames(otu_mat)))
NROW(unique(colnames(otu_mat)))

############################################################################################################


# create a sample table

df_dummy <- df1
df_dummy$sample_name <- paste0(df_dummy$pig,"_",df_dummy$date,"_")
df_dummy <- df_dummy %>%
  dplyr::select(sample_name, everything()) 
df_dummy <- df_dummy[,1:2]
head(df_dummy)

df_dummy$sample_name2 <- df_dummy$sample_name
df_dummy <- cSplit(df_dummy, "sample_name2", sep="_")

colnames(df_dummy) <- c("sample_name","cohort","pig","date")
NROW(df_dummy)

df_dummy$sample_name <- df_dummy$sample_name %<>%
  gsub('_$', '', .)  # removes the last _

df_dummy <- unique(df_dummy)
sample_df <- as.data.frame(df_dummy)

rownames(sample_df) <- sample_df[,1]

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
                                        "Dscour",
                                        "ColiGuard", 
                                        "Neomycin",
                                        "NeoD",
                                        "NeoC"))

############################################################################################################


# create phyloseq object

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)

carbom <- phyloseq(OTU,TAX,samples)
carbom


############################################################################################################

# necessary to rarefy? if yes, I ll need to convert counts to integers first 
#rarecurve(t(otu_table(carbom)), step=50, cex=0.5)

# SUBSETTING phyloseq obejct

# Keep only samples to be analyzed
#carbom <- subset_samples(carbom, date =="t2")
carbom <- subset_samples(carbom, (date %in% c("t0","t2","t4","t6","t8","t10")))

#subset to what you want
unique(taxa_mat_df$phylum)
carbom <- subset_taxa(carbom, (phylum %in% c("p__Firmicutes","p__Tenericutes","p__Actinobacteria",
                                             "p__Proteobacteria","p__Bacteroidetes","p__Spirochaetes",
                                             "p__Euryarchaeota","p__Chlamydiae","p__Synergistetes",
                                             "p__Verrucomicrobia")))


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

# BAR PLOT

pdf("phyloseq_abundance_barplot.pdf")
# BAR GRAPH - all samples 
plot_bar(carbom, fill = "phylum") + 
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
  theme(axis.text.x = element_blank())
# BAR GRAPH - by time point
plot_bar(carbom, fill = "phylum") + 
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
# BAR GRAPH - by time point - class
plot_bar(carbom, fill = "class") + 
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
dev.off()

######################

# HEATMAP

# plot_heatmap(carbom, method = "NMDS", distance = "bray")

# keep only very abundant OTUs
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.20) > 0, TRUE)

# HEATMAP with only most abundant OTUs
# plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

# HEATMAP with only most abundant OTUs - with names 
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "class", taxa.order = "class", 
             trans=NULL, low="beige", high="red", na.value="beige")

######################

# DIVERSITY 

pdf("phyloseq_diversity.pdf")
plot_richness(carbom, measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), x="cohort", color="date")
dev.off()


######################

# ORDINATION 

carbom.ord <- ordinate(carbom, "NMDS", "bray")

pdf("phyloseq_ordination.pdf")
plot_ordination(carbom, carbom.ord, type="samples", color="date", #shape= "cohort", 
                title="OTUs") + 
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

pdf("phyloseq_network.pdf")
ig = make_network(carbom_abund, type = "samples", distance = "bray", max.dist = 0.3)
plot_network(ig, carbom_abund, color = "date", shape = "cohort", line_weight = 0.3, 
             label = NULL, title = "sample network - Bray-Curtis distance")
dev.off()

############################################################################################################

# DESEQ2 (to be continued)

library(microbiome)
library(DESeq2)

d <- carbom

otu_table(d)
tax_table(d)
sample_data(d)

# Start by converting phyloseq object to deseq2 format
ds2 <- phyloseq_to_deseq2(d, ~ date + cohort)

# Run DESeq2 analysis (all taxa at once!)
dds <- DESeq(ds2)

# Investigate results
res <- results(dds)
deseq.results <- as.data.frame(res)
df <- deseq.results
df$taxon <- rownames(df)
df <- df %>% arrange(log2FoldChange, padj)

# Print the results; flitered and sorted by pvalue and effectsize
library(knitr)
df <- df %>% filter(pvalue < 0.05 & log2FoldChange > 1.5) %>%
  arrange(pvalue, log2FoldChange)
kable(df, digits = 5)


# For comparison purposes, assess significances and effect sizes based on Wilcoxon test.
# https://microbiome.github.io/tutorials/all.html


############################################################################################################



