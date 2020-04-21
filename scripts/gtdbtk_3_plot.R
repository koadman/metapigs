
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
        title="Phyla distribution from all MAGs (gtdbtk) - most abundant", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# least abundant 
treemap(phylum_counts_least_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "perc",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (gtdbtk) - least abundant", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()


