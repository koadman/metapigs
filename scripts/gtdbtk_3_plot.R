
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
        title="Phyla distribution from all MAGs (GTDB) - most common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# least abundant 
treemap(phylum_counts_least_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "perc",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (GTDB) - least common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()




######################################################################



# Species overall - based on gtdbtk assignment



species_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(family,species)]

# all 
species_counts <- species_counts %>%
  #filter(!Freq<10) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  arrange(desc(perc))
species_counts$label <- paste(paste(species_counts$species,
                                    species_counts$perc),"%")

# most abundant 
species_counts_most_ab <- species_counts[1:50]


pdf("gt_treemap_species.pdf")
# most abundant 
treemap(species_counts_most_ab, #Your data frame object
        index=c("family","label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Species distribution from all MAGs (GTDB) - 50 most common", #Customize your title
        fontsize.title = 15, #Change the font size of the title
        overlap.labels=0,
        fontsize.labels=c(13,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","black"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("left", "top"), 
          c("right", "bottom")
        )
        #fontsize.labels = 8
)
dev.off()


######################################################################



# Genus overall - based on gtdbtk assignment



genus_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(order,genus)]

# all 
genus_counts <- genus_counts %>%
  #filter(!Freq<10) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  arrange(desc(perc))
genus_counts$label <- paste(paste(genus_counts$genus,
                                  genus_counts$perc),"%")

# most abundant 
genus_counts_most_ab <- genus_counts[1:50]


pdf("gt_treemap_genus.pdf")
# most abundant 
treemap(genus_counts_most_ab, #Your data frame object
        index=c("order","label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Genus distribution from all MAGs (GTDB) - 50 most common", #Customize your title
        fontsize.title = 15, #Change the font size of the title
        overlap.labels=0,
        fontsize.labels=c(13,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","black"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("left", "top"), 
          c("right", "bottom")
        )
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


# SPECIES

# STEP 2.

# sum all the norm values that fall within same pig,date,taxa_2
df3 <- df2 %>%
  dplyr::group_by(pig,date,species) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample =paste0(df3$date,"_",df3$pig)
df3$date <- NULL
df3$pig <- NULL

# long to wide format
df4 <- df3 %>%
  pivot_wider(names_from = species, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 
head(df4)


#################################



# get a quick cohorts to pig table
cohorts <- df %>% dplyr::select(cohort,pig,date) %>% distinct()
cohorts$sample <- paste0(cohorts$date,"_",cohorts$pig)
cohorts <- as.data.frame(cohorts)


df5 <- inner_join(cohorts,df4) 
df5$sample <- paste0(df5$date,"_",df5$cohort)

df5$pig <- NULL
df5$date <- NULL
df5$cohort <- NULL



df6 <- df5 %>%
  group_by(sample) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)


df6 <- as.data.frame(df6)
rowSums(df6[,-1])



rownames(df6) <- df6$sample
df6$sample <- NULL
m <- as.matrix(df6)

df6.pca <- prcomp(m, center = FALSE,scale. = FALSE)
summary(df6.pca)

# to get samples info showing on PCA plot
this_mat_samples <- data.frame(sample=rownames(m)) 
this_mat_samples <- cSplit(indt = this_mat_samples, "sample", sep = "_", drop = NA)

# reorder dates 
this_mat_samples$sample_1  = factor(this_mat_samples$sample_1, levels=c("t0",
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

gt_PC12 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample_2,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (1:2)) +
  theme_bw() +
  xlim(c(-2,1)) +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))
gt_PC34 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample_2,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (3:4)) +
  theme_bw() +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))


gt_PCA <- ggarrange(gt_PC12,gt_PC34,
                    ncol=2,legend = "right",
                    common.legend=TRUE)

pdf("gt_PCA.pdf", width=7,height=4)
gt_PCA
dev.off()
