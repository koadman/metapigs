
library(readr)
library(splitstackshape)
library(treemap)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(robCompositions)
library(ggbiplot)

setwd("~/Desktop/metapigs_dry/gtdbtk")
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


######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv("gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################

# taxa overall - based on gtdbtk assignment


phylum_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(phylum)]

# most abundant 
phylum_counts_most_ab <- phylum_counts %>%
  filter(!Freq<2) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  # if percentage lower than 1.4 group as "others"
  mutate(phylum = ifelse(perc < 0.7, "others", phylum)) %>%
  group_by(phylum) %>%
  dplyr::summarise(perc=sum(perc))
  
# least abundant
phylum_counts_least_ab <- phylum_counts %>%
  filter(!Freq<2) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  filter(perc<0.7)

phylum_counts_most_ab$label <- paste(paste(phylum_counts_most_ab$phylum,
                                           phylum_counts_most_ab$perc,sep = "\n"),"%")
phylum_counts_least_ab$label <- paste(paste(phylum_counts_least_ab$phylum,
                                            phylum_counts_least_ab$perc,sep = "\n"),"%")


pdf("gt_treemap_phyla.pdf")
# most abundant 
treemap(phylum_counts_most_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "perc",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (gtdbtk) - most common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# least abundant 
treemap(phylum_counts_least_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "perc",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (gtdbtk) - least common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()




######################################################################



# Species overall - based on gtdbtk assignment


species_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(species)]

# all 
species_counts <- species_counts %>%
  #filter(!Freq<10) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  arrange(desc(perc))
species_counts$label <- paste(paste(species_counts$species,
                                    species_counts$perc,sep = "\n"),"%")

# most abundant 
species_counts_most_ab <- species_counts[1:50]

# least abundant
species_counts_least_ab <- species_counts[51:100]


pdf("gt_treemap_species.pdf")
# all 
treemap(species_counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Species distribution from all MAGs (gtdbtk) - all", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# most abundant 
treemap(species_counts_most_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Species distribution from all MAGs (gtdbtk) - 50 most common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# least abundant 
treemap(species_counts_least_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Species distribution from all MAGs (gtdbtk) - from 50th to 100th most common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()




######################################################################



# Genus overall - based on gtdbtk assignment


genus_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(genus)]

# all 
genus_counts <- genus_counts %>%
  #filter(!Freq<10) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  arrange(desc(perc))
genus_counts$label <- paste(paste(genus_counts$genus,
                                  genus_counts$perc,sep = "\n"),"%")

# most abundant 
genus_counts_most_ab <- genus_counts[1:50]

# least abundant
genus_counts_least_ab <- genus_counts[51:100]


pdf("gt_treemap_genus.pdf")
# all 
treemap(genus_counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="genus distribution from all MAGs (gtdbtk) - all", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# most abundant 
treemap(genus_counts_most_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="genus distribution from all MAGs (gtdbtk) - 50 most common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# least abundant 
treemap(genus_counts_least_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="genus distribution from all MAGs (gtdbtk) - from 50th to 100th most common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()



######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################



# Principal component analysis: clustering of bins based on phyla with time (labels per cohort)




# merge info 

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


df <- df0 %>%
  filter(!cohort=="Mothers")


#################################
# STEP 1.

# normalization for library size 
df2 <- df %>%
  dplyr::group_by(pig,date) %>%
  dplyr::mutate(norm_value=value/sum(value)) 
NROW(df2)
head(df2)

# # test:
# test <- df2 %>%
#   filter(pig=="14159") %>%
#   filter(date=="t0") %>%
#   dplyr::mutate(norm_value=value/sum(value))
# head(test)
# sum(test$norm_value)

#################################

# get a quick cohorts to pig table to use later for labelling
cohorts <- df0 %>% dplyr::select(cohort,pig) %>% distinct()

df2 
# ready to use

######################################################################################################
######################################################################################################


# PHYLUM

# STEP 2.

# sum all the norm values that fall within same pig,date,taxa_2
df3 <- df2 %>%
  dplyr::group_by(pig,date,phylum) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample =paste0(df3$date,"_",df3$pig)
df3$date <- NULL
df3$pig <- NULL

# long to wide format
df4 <- df3 %>%
  pivot_wider(names_from = phylum, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 

sample <- data.frame(sample = df4[,1])
df4 <- df4[,-1]

my_minimum <- (min(df4[df4 > 0])*0.1)*0.1

df4 <- df4+my_minimum

rownames(df4) <- sample$sample

#################

rowSums(df4)
df4_eclr <- cenLR(df4)
clr_norm_df <- df4_eclr$x.clr

# run PCA
df4.pca <- prcomp(clr_norm_df, center = TRUE,scale. = FALSE)

PCA_phylum <- ggbiplot(df4.pca,
                        var.axes = TRUE, labels = NULL,
                        groups=substr(rownames(clr_norm_df),1,3),
                        ellipse=TRUE,
                        choices = (1:2)) +
  theme_minimal() +
  ggtitle("Principal component analysis (cenLR) based on phylum (gtdbtk)")

###################################################################################################

# CLASS

df3 <- df2 %>%
  dplyr::group_by(pig,date,class) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample =paste0(df3$date,"_",df3$pig)
df3$date <- NULL
df3$pig <- NULL

df4 <- df3 %>%
  pivot_wider(names_from = class, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 

sample <- data.frame(sample = df4[,1])
df4 <- df4[,-1]
my_minimum <- (min(df4[df4 > 0])*0.1)*0.1
df4 <- df4+my_minimum
rownames(df4) <- sample$sample

rowSums(df4)
df4_eclr <- cenLR(df4)
clr_norm_df <- df4_eclr$x.clr

# run PCA
df4.pca <- prcomp(clr_norm_df, center = TRUE,scale. = FALSE)

PCA_class <- ggbiplot(df4.pca,
                      var.axes = TRUE, labels = NULL,
                      groups=substr(rownames(clr_norm_df),1,3),
                      ellipse=TRUE,
                      choices = (1:2)) +
  theme_minimal() +
  ggtitle("Principal component analysis (cenLR) based on class (gtdbtk)")

######################################################################################################

# ORDER

df3 <- df2 %>%
  dplyr::group_by(pig,date,order) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample =paste0(df3$date,"_",df3$pig)
df3$date <- NULL
df3$pig <- NULL

df4 <- df3 %>%
  pivot_wider(names_from = order, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 

sample <- data.frame(sample = df4[,1])
df4 <- df4[,-1]
my_minimum <- (min(df4[df4 > 0])*0.1)*0.1
df4 <- df4+my_minimum
rownames(df4) <- sample$sample

rowSums(df4)
df4_eclr <- cenLR(df4)
clr_norm_df <- df4_eclr$x.clr

# run PCA
df4.pca <- prcomp(clr_norm_df, center = TRUE,scale. = FALSE)

PCA_order <- ggbiplot(df4.pca,
                      var.axes = TRUE, labels = NULL,
                      groups=substr(rownames(clr_norm_df),1,3),
                      ellipse=TRUE,
                      choices = (1:2)) +
  theme_minimal() +
  ggtitle("Principal component analysis (cenLR) based on order (gtdbtk)")



######################################################################################################

# FAMILY

df3 <- df2 %>%
  dplyr::group_by(pig,date,family) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample =paste0(df3$date,"_",df3$pig)
df3$date <- NULL
df3$pig <- NULL

df4 <- df3 %>%
  pivot_wider(names_from = family, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 

sample <- data.frame(sample = df4[,1])
df4 <- df4[,-1]
my_minimum <- (min(df4[df4 > 0])*0.1)*0.1
df4 <- df4+my_minimum
rownames(df4) <- sample$sample

rowSums(df4)
df4_eclr <- cenLR(df4)
clr_norm_df <- df4_eclr$x.clr

# run PCA
df4.pca <- prcomp(clr_norm_df, center = TRUE,scale. = FALSE)

PCA_family <- ggbiplot(df4.pca,
                       var.axes = FALSE, labels = NULL,
                       groups=substr(rownames(clr_norm_df),1,3),
                       ellipse=TRUE,
                       choices = (1:2)) +
  theme_minimal() +
  ggtitle("Principal component analysis (cenLR) based on family (gtdbtk)")


######################################################################################################

# GENUS

df3 <- df2 %>%
  dplyr::group_by(pig,date,genus) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample =paste0(df3$date,"_",df3$pig)
df3$date <- NULL
df3$pig <- NULL

df4 <- df3 %>%
  pivot_wider(names_from = genus, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 

sample <- data.frame(sample = df4[,1])
df4 <- df4[,-1]
my_minimum <- (min(df4[df4 > 0])*0.1)*0.1
df4 <- df4+my_minimum
rownames(df4) <- sample$sample

rowSums(df4)
df4_eclr <- cenLR(df4)
clr_norm_df <- df4_eclr$x.clr

# run PCA
df4.pca <- prcomp(clr_norm_df, center = TRUE,scale. = FALSE)

PCA_genus <- ggbiplot(df4.pca,
                      var.axes = FALSE, labels = NULL,
                      groups=substr(rownames(clr_norm_df),1,3),
                      ellipse=TRUE,
                      choices = (1:2)) +
  theme_minimal() +
  ggtitle("Principal component analysis (cenLR) based on genus (gtdbtk)")



######################################################################################################

# species

df3 <- df2 %>%
  dplyr::group_by(pig,date,species) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample =paste0(df3$date,"_",df3$pig)
df3$date <- NULL
df3$pig <- NULL

df4 <- df3 %>%
  pivot_wider(names_from = species, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 

sample <- data.frame(sample = df4[,1])
df4 <- df4[,-1]
my_minimum <- (min(df4[df4 > 0])*0.1)*0.1
df4 <- df4+my_minimum
rownames(df4) <- sample$sample

rowSums(df4)
df4_eclr <- cenLR(df4)
clr_norm_df <- df4_eclr$x.clr

# run PCA
df4.pca <- prcomp(clr_norm_df, center = TRUE,scale. = FALSE)

PCA_species <- ggbiplot(df4.pca,
                        var.axes = FALSE, labels = NULL,
                        groups=substr(rownames(clr_norm_df),1,3),
                        ellipse=TRUE,
                        choices = (1:2)) +
  theme_minimal() +
  ggtitle("Principal component analysis (cenLR) based on species (gtdbtk)")


pdf("gt_PCA.pdf")
PCA_phylum
PCA_class
PCA_order
PCA_family
PCA_genus
PCA_species
dev.off()


