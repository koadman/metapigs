
library(readr)
library(splitstackshape)
library(treemap)
library(dplyr)
library(data.table)

setwd("~/Desktop/metapigs_dry/gtdbtk")
basedir = "~/Desktop/metapigs_dry/gtdbtk/"



# load gtdbtk assignments of the bins
gtdbtk_bins <- read_delim(paste0(basedir,"all_concatenated_220_essential.kraken"), 
                                             "\t", escape_double = FALSE, col_names = FALSE, 
                                             trim_ws = TRUE)

gtdbtk_bins <- gtdbtk_bins %>% 
  mutate(X3 = gsub("taxid.*", "", X3)) %>% 
  mutate(X2 = gsub(".fa", "", X2)) 

gtdbtk_bins <- cSplit(gtdbtk_bins, "X3"," (")
gtdbtk_bins$X1 <- NULL
gtdbtk_bins$X4 <- NULL

gtdbtk_bins <- cSplit(gtdbtk_bins, "X2","_")
colnames(gtdbtk_bins) <- c("species","pig","bin")



# load gtdb dictionary
bac_arc <- read_csv(paste0(basedir,"bac120_arc122_dictionary"),col_names = TRUE)


# SPECIES

# retrieve upper level taxa of species that have been assigned to my bins
named_byspecies <- inner_join(gtdbtk_bins,bac_arc)
NROW(named_byspecies) # were succesfully found as species (out of 51175)
NROW(unique(named_byspecies$species)) # <- these are the numbers of unique bins that have been assigned upper taxonomic levels
NROW(unique(gtdbtk_bins$species)) # <- these are the numbers of unique bins 

NOTnamed_byspecies <- anti_join(gtdbtk_bins,bac_arc)
colnames(NOTnamed_byspecies) <- c("genus","pig","bin")
NROW(NOTnamed_byspecies) # total number of bins (out of 51175) not found as species as these were not classified as species to start with
View(NOTnamed_byspecies)

# GENUS
bac_arcG <- bac_arc %>% 
  group_by(genus) %>% 
  filter(row_number()==1)
named_bygenus <- inner_join(NOTnamed_byspecies,bac_arcG)
named_bygenus$species <- NA
NROW(named_bygenus)
NROW(unique(named_bygenus$genus)) 

# FAMILY
bac_arcF <- bac_arc %>% 
  group_by(family) %>% 
  filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("family","pig","bin")
named_byfamily <- inner_join(NOTnamed_byspecies,bac_arcF)
named_byfamily$species <- NA
named_byfamily$genus <- NA
NROW(named_byfamily)
NROW(unique(named_byfamily$family)) 


# ORDER
bac_arcO <- bac_arc %>% 
  group_by(order) %>% 
  filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("order","pig","bin")
named_byorder <- inner_join(NOTnamed_byspecies,bac_arcO)
named_byorder$species <- NA
named_byorder$genus <- NA
named_byorder$family <- NA
NROW(named_byorder)
NROW(unique(named_byorder$order)) 

# CLASS
bac_arcC <- bac_arc %>% 
  group_by(class) %>% 
  filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("class","pig","bin")
named_byclass <- inner_join(NOTnamed_byspecies,bac_arcC)
named_byclass$species <- NA
named_byclass$genus <- NA
named_byclass$family <- NA
named_byclass$order <- NA
NROW(named_byclass)
NROW(unique(named_byclass$class)) 


# PHYLUM
bac_arcP <- bac_arc %>% 
  group_by(phylum) %>% 
  filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("phylum","pig","bin")
named_byphylum <- inner_join(NOTnamed_byspecies,bac_arcP)
named_byphylum$species <- NA
named_byphylum$genus <- NA
named_byphylum$family <- NA
named_byphylum$order <- NA
named_byphylum$class <- NA
NROW(named_byphylum)
NROW(unique(named_byphylum$phylum)) 

# domain
bac_arcD <- bac_arc %>% 
  group_by(domain) %>% 
  filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("domain","pig","bin")
named_bydomain <- right_join(NOTnamed_byspecies,bac_arcD)
named_bydomain$species <- NA
named_bydomain$genus <- NA
named_bydomain$family <- NA
named_bydomain$order <- NA
named_bydomain$class <- NA
named_bydomain$phylum <- NA
NROW(named_bydomain)
NROW(unique(named_bydomain$domain)) 

all <- rbind(named_byspecies,
      named_bygenus,
      named_byfamily,
      named_byorder,
      named_byclass,
      named_byphylum,
      named_bydomain)
NROW(all)

NROW(gtdbtk_bins)-NROW(all)

NROW(unique(paste0(gtdbtk_bins$pig,gtdbtk_bins$bin)))
NROW(unique(paste0(all$pig,all$bin)))

View(all)


fwrite(x = all, file = "gtdbtk_bins_completeTaxa")



unique(all$phylum)

# this is how many bins have been classified at phylum level 
NROW(which(!is.na(all$phylum)))/NROW(all)*100

# this is how many bins have been classified at class level 
NROW(which(!is.na(all$class)))/NROW(all)*100

# this is how many bins have been classified at order level 
NROW(which(!is.na(all$order)))/NROW(all)*100

# this is how many bins have been classified at family level 
NROW(which(!is.na(all$family)))/NROW(all)*100

# this is how many bins have been classified at genus level 
NROW(which(!is.na(all$genus)))/NROW(all)*100

# this is how many bins have been classified at species level 
NROW(which(!is.na(all$species)))/NROW(all)*100






phylum_counts <- setDT(all)[, .(Freq = .N), by = .(phylum)]

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


pdf("treemap_gtdbtk_phyla_new.pdf")
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



phylum_counts <- setDT(all)[, .(Freq = .N), by = .(order)]
