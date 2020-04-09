
library(readr)
library(splitstackshape)
library(treemap)
library(dplyr)
library(data.table)


setwd("~/Desktop/metapigs_dry/gtdbtk")
basedir = "~/Desktop/metapigs_dry/gtdbtk/"

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv("gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


unique(gtdbtk_bins$phylum)

# this is how many bins have been classified at phylum level 
NROW(which(!is.na(gtdbtk_bins$phylum)))/NROW(gtdbtk_bins)*100

# this is how many bins have been classified at class level 
NROW(which(!is.na(gtdbtk_bins$class)))/NROW(gtdbtk_bins)*100

# this is how many bins have been classified at order level 
NROW(which(!is.na(gtdbtk_bins$order)))/NROW(gtdbtk_bins)*100

# this is how many bins have been classified at family level 
NROW(which(!is.na(gtdbtk_bins$family)))/NROW(gtdbtk_bins)*100

# this is how many bins have been classified at genus level 
NROW(which(!is.na(gtdbtk_bins$genus)))/NROW(gtdbtk_bins)*100

# this is how many bins have been classified at species level 
NROW(which(!is.na(gtdbtk_bins$species)))/NROW(gtdbtk_bins)*100






phylum_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(phylum)]

phylum_counts_most_ab <- phylum_counts %>%
  filter(!Freq<2) %>%
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  # most abundant 
  filter(perc>1.4)
phylum_counts_least_ab <- phylum_counts %>%
  filter(!Freq<2) %>%
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  # least abundant 
  filter(perc<1.4)

phylum_counts_most_ab$label <- paste(paste(phylum_counts_most_ab$phylum,
                                           phylum_counts_most_ab$perc,sep = "\n"),"%")
phylum_counts_least_ab$label <- paste(paste(phylum_counts_least_ab$phylum,
                                            phylum_counts_least_ab$perc,sep = "\n"),"%")


pdf("gt_treemap_phyla.pdf")
# most abundant 
treemap(phylum_counts_most_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (gtdbtk db) - most abundant", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# least abundant 
treemap(phylum_counts_least_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (gtdbtk db) - least abundant", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()


