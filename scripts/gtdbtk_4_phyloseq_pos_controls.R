
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
# no_reps_pos_controls.csv (BINS COUNTS)


# OUTPUTS:
# gt_phylo_PosControls_barplot.pdf
# gt_phylo_PosControls_heatmap.pdf



######################################################################

# counts data 

no_reps_pos_controls <- read.csv(paste0(basedir,"pos_controls/no_reps_pos_controls.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension to match bins in checkm df 
no_reps_pos_controls$bin <- gsub(".fa","", no_reps_pos_controls$bin)
head(no_reps_pos_controls)
NROW(no_reps_pos_controls)

colnames(no_reps_pos_controls)[colnames(no_reps_pos_controls) == 'isolation_source'] <- 'pig'

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
NROW(no_reps_pos_controls)
head(gtdbtk_bins)
head(no_reps_pos_controls)
df <- merge(no_reps_pos_controls, gtdbtk_bins, by=c("pig","bin"))

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
keep <- c("cohort","pig","bin","variable","value","gOTU")
df1 <- df[ , (names(df) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$variable)

NROW(df1)
# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,variable,gOTU) %>%
  dplyr::summarise(all_bins_value = sum(value))

NROW(df2)
NROW(unique(paste0(df2$pig,df2$variable)))


# assign a unique sample name 
df2$sample <- paste0(df2$pig,"_",df2$variable)
# remove now pig and date (redundant)
df2$pig <- NULL
df2$variable <- NULL

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

sample_df$sample <- paste0(sample_df$pig,"_",sample_df$variable)
NROW(unique(sample_df$sample))

sample_df <- sample_df %>%
  dplyr::select(sample,pig,variable,cohort) %>%
  group_by(sample) %>%
  slice(1)

sample_df$gOTU <- NULL
sample_df <- as.data.frame(sample_df)

NROW(sample_df)
head(sample_df)

rownames(sample_df) <- sample_df[,1]
# ready


######################################################################


# create phyloseq object

OTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)

carbom <- phyloseq(OTU,TAX,samples)


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
sampleOrder = sort(sample_names(carbom))


# HEATMAP time - genus, family, order, etc ...
pdf("gt_phylo_PosControls_heatmap.pdf")
plot_heatmap(carbom, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  ggtitle(label = "Positive Controls Composition (gtdbtk)") 
dev.off()

######################


# BAR PLOT

# BAR GRAPH - by time point
pdf("gt_phylo_PosControls_barplot.pdf")
plot_bar(carbom, fill = "species") +
  geom_bar(aes(color=species), stat="identity", position="stack") +
  facet_grid(~pig,scales="free_x") +
  theme(axis.text.x = element_text()) +
  ggtitle(label = "Positive Controls Composition (gtdbtk)") 
dev.off()


######################################################################
######################################################################
