

# upload libraries
library(tidyverse)
library(gplots)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(splitstackshape)
library(treemap)
library(dplyr)
library(data.table)
library(readxl)
library(pheatmap)
library(robCompositions)
library(ggbiplot)

basedir = "/Users/12705859/Desktop/metapigs_dry/"
setwd("/Users/12705859/Desktop/metapigs_dry/checkm/")

# input files: 
# all_checkm_output.tsv 
# checkm/checkm_all_nearly
# no_reps_all.csv (BINS COUNTS)


# OUTPUTS:


# about bins in general :
#  - numbers.txt

# from all_checkm_output :
#  - cm_Compl_vs_contam.pdf
#  - cm_numbers.txt
#  - cm_contigs_distribution.pdf
#  - cm_scaffolds_predicted_genes.pdf

# from BINS COUNTS + checkm_all_nearly: 
# based on bins frequency: 
## - cm_treemap_phylum.pdf
## - cm_treemap_class.pdf
## - cm_treemap_order.pdf
## - numbers.txt
# based on (counts) log10 relative abundance: 
## - cm_parallel_coordinates_phyla.pdf
# based on (counts) relative abundance: 
## - cm_rel_ab_phyla.pdf
## - cm_rel_ab_phyla_cohorts.pdf
## - cm_balloonplot_phyla.pdf
## - cm_balloonplot_phyla_cohorts.pdf
# based on (counts) cenLR relative abundance: 
## - cm_PCA.pdf



######################################################################

# template to collect and store bins info 

sink("numbers.txt")
start_message <- " ########################## BINS NUMBERS ########################## "
start_message
sink()

# template to collect and store checkM info 

sink("cm_numbers.txt")
start_message <- " ########################## CHECKM ANALYSIS ########################## "
start_message
sink()

######################################################################

# checkM output of ALL bins 

# upload file
# careful cause we don't have the pigID to distinguish which bins to which sample
all_checkm_output <- read_delim(paste0(basedir,"checkm/all_checkm_output.tsv"),
                                "\t", escape_double = FALSE, trim_ws = TRUE)

all_checkm_output <- dplyr::filter(all_checkm_output, !grepl("Completeness",Completeness))
all_checkm_output$Completeness <- as.numeric(all_checkm_output$Completeness)
all_checkm_output$Contamination <- as.numeric(all_checkm_output$Contamination)
NROW(all_checkm_output)


######################################################################


# checkM output of NEARLY COMPLETE bins 

# upload file
# careful cause we don't have the pigID to distinguish which bins to which sample
checkm_all_nearly <- read_delim("checkm_all_nearly", 
                                "\t", escape_double = FALSE, col_types = cols(pigid = col_character()), 
                                trim_ws = TRUE)



######################################################################


# upload bins with counts (sample-dereplicated- output of 7.R)

no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)

head(no_reps_all)


######################################################################


# NUMBERS 


##########################

# Bins distribution:
# how many bins each subject has:
az <- no_reps_all %>%
  dplyr::select(bin, pig, cohort) %>%
  distinct() 

az <- az %>%
  mutate(group = ifelse(cohort == "Mothers", "mothers (n=42)","piglets (n=126)"))

az <- az %>%
  group_by(group,pig) %>%
  dplyr::summarise(`number of metagenomes obtained`= n()) 
tail(az)

pdf("#bins_subject.pdf")
ggplot(data=az, mapping=aes(x=group, y=`number of metagenomes obtained`)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1)
dev.off()


##########################

# Bins per cohort: 
az2 <- no_reps_all %>%
  select(cohort,bin,pig) %>%
  group_by(cohort,pig) %>%
  dplyr::summarise(`number of metagenomes obtained`= n()) 
tail(az2)

pdf("#bins_cohort.pdf")
ggplot(data=az2, mapping=aes(x=cohort, y=`number of metagenomes obtained`)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1)
dev.off()

##########################


# Get the numbers: 
az <- no_reps_all %>%
  dplyr::select(bin,pig)
bz <- unique(az)
cz <- data.frame(table(bz$pig))

sink(file = "numbers.txt", 
     append = TRUE, type = c("output"))
paste0("total of bins (excludes the pos and neg controls bins): ", sum(cz$Freq) )
sink()


##########################


# number of samples per time point: 

az <- no_reps_all %>%
  filter(!cohort=="Mothers") %>%
  dplyr::select(pig, date) %>%
  distinct() 

az <- az %>%
  group_by(date) %>%
  dplyr::summarise(`number of samples`= n()) 
tail(az)

# reorder dates 
az$date  = factor(az$date, levels=c("t0",
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
az <- az %>%
  mutate(sampling = ifelse(date == "t0" | date == "t2" | 
                             date == "t4" | date == "t6" | 
                             date == "t8" | date == "t10" , "all subjects", "subset"))

pdf("#samples_per_timepoint.pdf")
ggplot(data=az, mapping=aes(x=date, y=`number of samples`, color=sampling)) + 
  geom_point() +
  geom_line(aes(group = sampling), linetype = 2) +
  theme_bw() + 
  labs(title="Number of collected samples per timepoint", 
       x = "timepoint",
       y = "number of samples",
       subtitle=NULL)
dev.off()


##########################


# Correlation between timepoints available and bins obtained per subject: 
# How many timepoints scale to how many bins? 

# BINS PER SUBJ
az1 <- no_reps_all %>%
  dplyr::select(bin, pig) %>%
  distinct() %>%
  group_by(pig) %>%
  dplyr::summarise(`number of metagenomes per subject`= n()) 
tail(az1)
NROW(az1)
head(az1)

sink(file = "numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Bins per subject: ")
summary(az1$`number of metagenomes per subject`)
sink()

# TIMEPOINTS PER SUBJ
az2 <- no_reps_all %>%
  dplyr::select(date, pig) %>%
  distinct() %>%
  group_by(pig) %>%
  dplyr::summarise(`number of timepoints per subject`= n()) 
tail(az2)
NROW(az2)

sink(file = "numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Timepoints per subject: ")
summary(az2$`number of timepoints per subject`)
sink()


az3 <- cbind(az1,az2)
az3$pig <- NULL

pdf("#bins_vs_#timepoints.pdf")
ggplot(data=az3, mapping=aes(x=`number of timepoints per subject`, y=`number of metagenomes per subject`)) + 
  geom_point(color='blue', size=0.6) +
  #geom_smooth(method = "lm", se = FALSE)+
  labs(title="Bins obtained versus number of timepoints available from each subject", 
       x = "number of timepoints available per subject",
       y = "number of bins obtained per subject",
       subtitle=NULL)
dev.off()



######################################################################

# bins distribution - cohorts & piglets

test <- no_reps_all %>%
  filter(!cohort=="Mothers") %>%
  dplyr::select(pig,bin,cohort) %>%
  distinct() %>%
  group_by(cohort,pig) %>%
  tally()  %>%
  mutate(perc=round(n/sum(n)*100,2)) 

# DT <- data.table(test)
# 
# DT[, id := seq_len(.N), by = cohort]
# DT[, id := rowid(cohort)]
# 
# DT <- DT %>%
#   mutate(newID = paste0("s",id))

pdf("bins_distribution_over_cohorts_piglets.pdf")
treemap(test, 
        index=c("cohort","pig"), 
        vSize="perc", 
        type="index",                            # How you color the treemap. type help(treemap) for more info
        palette = "Set1",                        # Select your color palette from the RColorBrewer presets or make your own.
        fontsize.title=12,                       # Size of the title
        title=paste0("Bins distribution over cohorts and piglets (",NROW(unique(paste0(no_reps_all$pig,no_reps_all$bin))),")")
)
dev.off()


######################################################################


# preparing data for Completeness vs Contamination dot plot

# keep (for the plot) only bins with >70% completeness and <10% contamination
df_70 <- all_checkm_output %>%
  filter(Completeness >= 70)
df_70_10 <- df_70 %>%
  filter(Contamination <= 10)

# assign name for plot legend 
df_70_10$type <- ifelse(df_70_10$Completeness >=90 & df_70_10$Contamination <=5, "high", "low")

# removing data less than <80 Complete to plot neatly
df_80_10 <- df_70_10 %>%
  filter(Completeness >= 80)
NROW(checkm_all_nearly)
df_80_10$type <- gsub("high",
                      paste0("Nearly complete genomes: ",
                             NROW(checkm_all_nearly),
                             " (",
                             round(NROW(checkm_all_nearly)/NROW(all_checkm_output)*100,2),
                             "%)"),
                      df_80_10$type)

medium_quality_almost <- df_70_10 %>%
  filter(Completeness >= 80) %>%
  filter(Completeness <= 100) %>%
  filter(Contamination <= 10)
medium_quality <- NROW(medium_quality_almost)-NROW(checkm_all_nearly)

df_80_10$type <- gsub("low",
                      paste0("Medium quality genomes: ",
                             medium_quality,
                             " (",
                             round(medium_quality/NROW(all_checkm_output)*100,2),
                             "%)"),
                      df_80_10$type)


# Assessment of genome quality
pdf("cm_Compl_vs_Contam.pdf")
ggplot(df_80_10, aes(x=Completeness, y=Contamination)) + 
  geom_point(aes(color=factor(type)),size=0.1,shape=21)+
  ylab("Contamination (%)")+
  xlab("Completeness (%)")+
  xlim(75,100)+
  theme_minimal()+
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.title = element_blank())
dev.off()

# similar to https://www.nature.com/articles/s41564-017-0012-7#Sec2


# library(scales)
# scales::hue_pal()(3)

# nearly complete bins 
df_90 <- all_checkm_output %>%
  filter(Completeness >= 90)
df_90_5 <- df_90 %>%
  filter(Contamination <= 5)

# perfect bins
df_99 <- all_checkm_output %>%
  filter(Completeness > 99)
df_perfect <- df_99 %>%
  filter(Contamination < 0.1)

sink(file = "cm_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Total number of bins is:   ",
       NROW(all_checkm_output))
paste0("Of which: >70% completeness and <10% contamination:   ",
       NROW(df_70_10)," ", round(NROW(df_70_10)/NROW(all_checkm_output),4)*100,"%")
paste0("Of which: >80% completeness and <10% contamination: ",
       NROW(df_80_10)," ", round(NROW(df_80_10)/NROW(all_checkm_output),4)*100,"%")
paste0("Of which: >=90% completeness and <=5% contamination: ",
       NROW(df_90_5)," ", round(NROW(df_90_5)/NROW(all_checkm_output),4)*100,"%")
paste0("Of which: >99% completeness and <0.1% contamination: ",
       NROW(df_perfect)," ", round(NROW(df_perfect)/NROW(all_checkm_output),4)*100,"%")
sink()


######################################################################


# Contig number and N50 distribution 

# as numeric 
df_90_5$`N50 (scaffolds)` <- as.numeric(df_90_5$`N50 (scaffolds)`)
df_90_5$`# contigs` <- as.numeric(df_90_5$`# contigs`)


pdf("cm_contigs_distribution.pdf")
# Distribution of contig number across bins 
hist(df_90_5$`# contigs`, 
     main = "Distribution of # of contigs across nearly complete genomes",
     breaks=100,
     xlab = "contigs")
hist(log10(df_90_5$`# contigs`), 
     main = "Distribution of # of contigs across nearly complete genomes",
     breaks=100,
     xlab = "log10(contigs)")
hist(df_90_5$`N50 (scaffolds)`, 
     main = "Distribution of N50 (scaffolds) across nearly complete genomes",
     breaks=100,
     xlab = "log10(N50 scaffolds)",
     xlim=c(0,3e+05))
hist(log10(df_90_5$`N50 (scaffolds)`), 
     main = "Distribution of N50 (scaffolds) across nearly complete genomes",
     breaks=100,
     xlab = "log10(N50 scaffolds)",
     xlim=c(3.5,6))
dev.off()


new <- all_checkm_output %>%
  select(`Bin Id`,`# genomes`,Completeness,Contamination, `# scaffolds`,`# predicted genes`)
new$id <- paste0(new$`Bin Id`,"_",new$`# genomes`)


new$`# scaffolds` <- as.numeric(new$`# scaffolds`)
new$`# predicted genes` <- as.numeric(new$`# predicted genes`)


new1 <- new %>%
  filter(Completeness >= 50) %>%
  filter(Completeness < 70) %>%
  filter(Contamination <= 10) %>%
  mutate(type="Compl50-70_Contam<10")

new2 <- new %>%
  filter(Completeness >= 70) %>%
  filter(Completeness < 90) %>%
  filter(Contamination <= 10) %>%
  mutate(type="Compl70-90_Contam<10")

new3 <- new %>%
  filter(Completeness >= 90) %>%
  filter(Contamination <= 10) %>%
  filter(Contamination > 5) %>%
  mutate(type="Compl>90_Contam5-10")

new4 <- new %>%
  filter(Completeness >= 90) %>%
  filter(Contamination <= 5) %>%
  mutate(type="Compl>90_Contam<5")

new0 <- new %>%
  filter(Completeness > 50) %>%
  filter(Contamination > 10) %>%
  mutate(type="Compl>50_Contam>10")

new00 <- new %>%
  filter(Completeness < 50) %>%
  mutate(type="Compl<50")

new_tot <- rbind(new00,new0,new1,new2,new3,new4)
NROW(new_tot)
NROW(all_checkm_output)

new23 <- rbind(new2,new3)

new23$type <- "Medium quality"

new1$type <- "Partial"

new4$type <- "Nearly complete"

new_plot <- rbind(new1,new23,new4)

new_plot$type  = factor(new_plot$type, levels=c("Nearly complete",
                                                "Medium quality",
                                                "Partial"))

pdf("cm_scaffolds_predicted_genes.pdf")
ggplot(data= new_plot, aes(x=`# scaffolds`, fill=type)) +
  geom_histogram(alpha=0.6, position = 'stack', binwidth=15) +
  scale_fill_manual(values=c("#45B4B8", #blue
                             "#F8766D", #red
                             "#6C6C6C")) + #grey
  theme_bw() +
  labs(fill="") +
  lims(x=c(0,750))
ggplot(data= new_plot, aes(x=`# predicted genes`, fill=type)) +
  geom_histogram(alpha=0.6, position = 'stack', binwidth=15) +
  scale_fill_manual(values=c("#45B4B8", #blue
                             "#F8766D", #red
                             "#6C6C6C")) + #grey
  theme_bw() +
  labs(fill="") +
  lims(x=c(0,5000))
dev.off()



######################################################################

# Nearly complete bins (>=90 <=5) taxonomic categorization: proportions explained at each taxonomic level: 

# -- proportion of >=90 and <=5 bins assigned to phylum/order/etc... : 

# how many of >90 <5 bins are classified at the each level? 
nn <- checkm_all_nearly %>% dplyr::rename(taxa = `Taxonomy (contained)`)
NROW(nn)

archaea <- filter(nn, grepl('k__Archaea', taxa))
bacteria <- filter(nn, grepl('k__Bacteria', taxa))
kingdom <- filter(nn, grepl('k__', taxa))
phylum <- filter(nn, grepl('p__', taxa))
class <- filter(nn, grepl('c__', taxa))
order <- filter(nn, grepl('o__', taxa))
family <- filter(nn, grepl('f__', taxa))
genus <- filter(nn, grepl('g__', taxa))
species <- filter(nn, grepl('s__', taxa))

sink(file = "cm_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("######################################################################")
paste0("Of the nearly complete bins: ")
paste0("Assigned to the Kingdom of Bacteria:   ",
       NROW(bacteria),"  ",round(NROW(bacteria)/NROW(nn),4)*100,"%")
paste0("Assigned to the Kingdom of Archaea:   ",
       NROW(archaea),"  ",round(NROW(archaea)/NROW(nn),4)*100,"%")
paste0("Resolved at Phylum level:   ",
       NROW(phylum),"  ",round(NROW(phylum)/NROW(nn),4)*100,"%")
paste0("Resolved at Class level:   ",
       NROW(class),"  ",round(NROW(class)/NROW(nn),4)*100,"%")
paste0("Resolved at Order level:   ",
       NROW(order),"  ",round(NROW(order)/NROW(nn),4)*100,"%")
paste0("Resolved at Family level:   ",
       NROW(family),"  ",round(NROW(family)/NROW(nn),4)*100,"%")
paste0("Resolved at Genus level:   ",
       NROW(genus),"  ",round(NROW(genus)/NROW(nn),4)*100,"%")
paste0("Resolved at Species level:   ",
       NROW(species),"  ",round(NROW(species)/NROW(nn),4)*100,"%")
sink()


######################################################################

# -- Proportions of ALL TAXA in the >90 <5 bins: 




k <- kingdom %>%
  mutate(taxa = sub(".*k__ *(.*?) *;p__.*", "\\1", taxa)) %>%
  mutate(taxa = gsub(".*k__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  mutate(perc=round(n/sum(n)*100,2)) %>%
  arrange(desc(perc))

p <- phylum %>%
  mutate(taxa = sub(".*p__ *(.*?) *;c__.*", "\\1", taxa)) %>%
  mutate(taxa = gsub(".*p__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  mutate(perc=round(n/sum(n)*100,2)) %>%
  arrange(desc(perc))

c <- class %>%
  mutate(taxa = sub(".*c__ *(.*?) *;o__.*", "\\1", taxa)) %>%
  mutate(taxa = gsub(".*c__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  mutate(perc=round(n/sum(n)*100,2)) %>%
  arrange(desc(perc))

o <- order %>%
  mutate(taxa = sub(".*o__ *(.*?) *;f__.*", "\\1", taxa)) %>%
  mutate(taxa = gsub(".*o__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  mutate(perc=round(n/sum(n)*100,2)) %>%
  arrange(desc(perc))

f <- family %>%
  mutate(taxa = sub(".*f__ *(.*?) *;g__.*", "\\1", taxa)) %>%
  mutate(taxa = gsub(".*f__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  mutate(perc=round(n/sum(n)*100,2)) %>%
  arrange(desc(perc))

g <- genus %>%
  mutate(taxa = sub(".*g__ *(.*?) *;s__.*", "\\1", taxa)) %>%
  mutate(taxa = gsub(".*g__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  mutate(perc=round(n/sum(n)*100,2)) %>%
  arrange(desc(perc))

s <- species %>%
  mutate(taxa = sub(".*s__ *(.*?) *", "\\1", taxa)) %>%
  mutate(taxa = gsub(".*s__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  mutate(perc=round(n/sum(n)*100,2)) %>%
  arrange(desc(perc))


p <- as.data.frame(p)
c <- as.data.frame(c)
o <- as.data.frame(o)
f <- as.data.frame(f)
g <- as.data.frame(g)
s <- as.data.frame(s)

sink(file = "cm_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("######################################################################")
paste0(" Kingdom distribution : ")
k
paste0("######################################################################")
paste0(" Phyla distribution : ")
p
paste0("######################################################################")
paste0(" Class distribution : ")
c
paste0("######################################################################")
paste0(" Order distribution : ")
o
paste0("######################################################################")
paste0(" Family distribution : ")
f
paste0("######################################################################")
paste0(" Genus distribution : ")
g
paste0("######################################################################")
paste0(" Species distribution : ")
s
sink()


######################################################################



# MERGE checkM info of nearly complete bins to bins counts : 

# rename cols of checkm nearly complete bins to match colnames of no_reps_all for dataframes to merge  
colnames(nn)[colnames(nn)=="pigid"] <- "pig"
colnames(nn)[colnames(nn)=="Bin Id"] <- "bin"
colnames(nn)[colnames(nn)=="Taxonomy (contained)"] <- "taxa"

nn <- nn %>%
  dplyr::select(pig,bin,taxa)

# merge 
df <- right_join(no_reps_all, nn, by=c("pig","bin"))
head(df)
NROW(df)


# split taxa column into several (kingdom, phylum, etc ...) 
df <- cSplit(df, "taxa", sep=";")

# reorder dates 
df$date  = factor(df$date, levels=c("tM", 
                                        "t0",
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
df$cohort  = factor(df$cohort, levels=c("Control", 
                                    "DScour",
                                    "ColiGuard", 
                                    "Neomycin",
                                    "NeoD",
                                    "NeoC"))



######################################################################


# treemap PHYLA (piglet samples, mothers exlcuded)

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_2) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_treemap[4] <- lapply(
  for_treemap[4], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

counts <- setDT(for_treemap)[, .(Freq = .N), by = .(taxa_2)]

counts <- counts %>%
  mutate(perc=round(Freq/sum(Freq)*100,2))

counts$label <- paste(paste(counts$taxa_2,counts$perc,sep = "\n"),"%")

pdf("cm_treemap_phyla.pdf")
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution - from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()

# treemap CLASS

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_3) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_treemap[4] <- lapply(
  for_treemap[4], 
  gsub, 
  pattern = "c__", 
  replacement = "", 
  fixed = TRUE)

counts <- setDT(for_treemap)[, .(Freq = .N), by = .(taxa_3)]

counts <- counts %>%
  mutate(perc=round(Freq/sum(Freq)*100,2))

counts$label <- paste(paste(counts$taxa_3,counts$perc,sep = "\n"),"%")

pdf("cm_treemap_class.pdf")
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Class distribution from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()

# treemap ORDER

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_4) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_treemap[4] <- lapply(
  for_treemap[4], 
  gsub, 
  pattern = "o__", 
  replacement = "", 
  fixed = TRUE)

counts <- setDT(for_treemap)[, .(Freq = .N), by = .(taxa_4)]

counts <- counts %>%
  mutate(perc=round(Freq/sum(Freq)*100,2))

counts$label <- paste(paste(counts$taxa_4,counts$perc,sep = "\n"),"%")

pdf("cm_treemap_order.pdf")
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Order distribution from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()




######################################################################################################

# parallel coordinates & relative abundance - PHYLA

# keep only rows that contain phyla info, discard others; also exclude mothers
for_parallel_coo <- df %>%
  dplyr::select(cohort,pig,bin,date,taxa_2,value) %>%
  filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_parallel_coo[5] <- lapply(
  for_parallel_coo[5], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

# general time change 
summs_for_parallel_coo <- for_parallel_coo %>% group_by(taxa_2,date) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 



pdf("cm_rel_ab_phyla.pdf")
ggplot(summs_for_parallel_coo, aes(fill=taxa_2, y=mean, x=date)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.title.y = element_text(),
        legend.title = element_text()) +
    labs(x = "collection date",
         y = "relative abundance (ratio)",
         fill = "Phylum")  
dev.off()


summs_for_parallel_coo_cohorts <- for_parallel_coo %>% group_by(taxa_2,date,cohort) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 

pdf("cm_rel_ab_phyla_cohorts.pdf")
ggplot(summs_for_parallel_coo_cohorts, aes(fill=taxa_2, y=mean, x=date)) + 
  geom_bar(position="fill", stat="identity")+
  facet_wrap(~cohort) +
  theme(axis.title.y = element_text(),
        legend.title = element_text()) +
  labs(x = "collection date",
       y = "relative abundance (ratio)",
       fill = "Phylum") 
dev.off()

# higher proportion of Tenericutes in Neomycin from t5 to t10 compared to other cohorts?
# Euryarchaeota "disapper" from t5 to t10. 
# Actinobacteria slowly increase in proportion with time
# igher bacteroidetes proportion between t3 and t4 in Neo than other cohorts? 
# Synergistetes prettu much disappear after t1 in Control, Dscour and ColiGuard,
# whereas they disappera later (after t2) in Neo and NeoD. 
# In NeoC in t4 still high proportion before gone. 
# Proteobacterai from t4 on higher proportion in Dscour, low in any other cohort. 


######################################################################################################

# PHYLA composition change through time: PARALLEL COORDINATES 

# keep only rows that contain phyla info, discard others; also exclude mothers
df1 <- df %>%
  dplyr::select(cohort,pig,bin,date,taxa_2,value) %>%
  filter(!cohort=="Mothers")

# remove "p__" before phylum 
df1[5] <- lapply(
  df1[5], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

df1 <- na.omit(df1)

df2 <- df1

# lib size normalization
df2 <- df2 %>% 
  group_by(pig,date) %>% 
  mutate(norm_value = (value/sum(value))) %>% 
  dplyr::select(pig, date, taxa_2, value, norm_value)

# sum all the norm values that fall within same pig,date,phylum:
df2 <- df2 %>%  
  group_by(pig,date,taxa_2) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))

# take the mean of each phylum by date:
df2 <- df2 %>%
  group_by(taxa_2,date) %>%
  dplyr::summarize(mean = mean(indiv_sum, na.rm=TRUE)) %>%
  mutate(perc=mean*100) %>%
  dplyr::select(taxa_2,date,perc)


df2 <- as.data.frame(df2)

pdf("cm_parallel_coordinates_phyla.pdf")
ggplot(df2, aes(x=date, y=log10(perc), group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "relative abundance (log10)",
       color = "Phylum")  +
  theme(legend.title = element_text()) 
dev.off()


######################################################################################################

# PHYLA composition change through time: BALLOONPLOTS - all


df3 <- df2
# pivot wider
df3 <- df3 %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

df4 <- df3
df4 <- as.data.frame(df4)
rownames(df4) <- df4[,1]
df4[,1] <- NULL
m <- as.table(as.matrix(df4))

pdf("cm_balloonplot_phyla.pdf")
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()


######################################################################################################

# PHYLA composition change through time: HEATMAP - all

pdf("cm_heatmap_phyla.pdf")
pheatmap(t(m), display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
dev.off()




######################################################################################################

# PHYLA composition change through time: BALLOONPLOTS - by cohort

df2 <- df1

# splitting into multiple dataframes (by cohort)
multi_df <- split( df2 , f = df2$cohort )

# construct an empty dataframe to build on 
final <- data.frame(
  pig = character(),
  bin = character(),
  date = character(),
  taxa_2 = character(),
  value = character(),
  stringsAsFactors = FALSE
)


for (single_df in multi_df) {
  
  single_df <- as.data.frame(single_df)
  coho <- as.character(single_df$cohort[1])
  
  # lib size normalization
  df2 <- single_df %>% 
    group_by(pig,date) %>% 
    mutate(norm_value = (value/sum(value))) %>% 
    dplyr::select(pig, date, taxa_2, value, norm_value)
  
  # sum all the norm values that fall within same pig,date,phylum:
  df2 <- df2 %>%  
    group_by(pig,date,taxa_2) %>%
    dplyr::summarise(indiv_sum = sum(norm_value))
  
  # take the mean of each phylum by date:
  df2 <- df2 %>%
    group_by(taxa_2,date) %>%
    dplyr::summarize(mean = mean(indiv_sum, na.rm=TRUE)) %>%
    mutate(perc=mean*100) %>%
    dplyr::select(taxa_2,date,perc)
  
  df2 <- as.data.frame(df2)
  df2$cohort <- coho
  
  final <- rbind(final,df2)
  
}

pdf("cm_balloonplot_phyla_cohorts.pdf")
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
# Control
fin <- final %>%
  dplyr::filter(cohort=="Control") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# DScour
fin <- final %>%
  dplyr::filter(cohort=="DScour") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# ColiGuard
fin <- final %>%
  dplyr::filter(cohort=="ColiGuard") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
# Neomycin
fin <- final %>%
  dplyr::filter(cohort=="Neomycin") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# NeoD
fin <- final %>%
  dplyr::filter(cohort=="NeoD") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# NeoC
fin <- final %>%
  dplyr::filter(cohort=="NeoC") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()



######################################################################################################
######################################################################################################

# Principal component analysis: clustering of bins based on phyla with time (labels per cohort)




# I hereby select taxa_2 only (corresponds to phylum) and remove any row where no phylum was resolved
df1 <- df %>%
  select(pig,bin,date,value,taxa_2,cohort)

df1 <- as.data.frame(na.omit(df1))

# remove "p__" before phylum
df1[5] <- lapply(
  df1[5],
  gsub,
  pattern = "p__",
  replacement = "",
  fixed = TRUE)

unique(df1$taxa_2)


#################################
# STEP 1.

# normalization for library size 
df2 <- df1 %>%
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
# STEP 2.

# sum all the norm values that fall within same pig,date,taxa_2
df3 <- df2 %>%
  dplyr::group_by(pig,date,taxa_2) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

# # test:
# test2 <- test %>%
#   group_by(taxa_2) %>%
#   dplyr::summarise(indiv_sum = sum(norm_value))
# head(test2)
# sum(test2$indiv_sum)

#################################
# STEP 3.

# long to wide format
df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 
head(df4)

# # test:
# test3 <- test2 %>%
#   pivot_wider(names_from = taxa_2, values_from = indiv_sum, values_fill = list(indiv_sum = 0))
# head(test3)
# sum(test3[1,])


#################################

# STEP 4. PCA, dots are cohort_date 

# get a quick cohorts to pig table 
cohorts <- df %>% dplyr::select(cohort,pig) %>% distinct()

# join the cohort info
df5 <- inner_join(df4,cohorts) %>%
  dplyr::mutate(coho_date_group=paste0(date,"_",cohort)) 
df5

# 
df6 <- df5 %>% 
  dplyr::group_by(coho_date_group) %>% 
  dplyr::summarise_if(is.numeric, funs(sum))
df6


rowSums(df6[,-1])
df6_eclr <- cenLR(df6[,-1])
clr_norm_df <- df6_eclr$x.clr

rownames(clr_norm_df) <- df6$coho_date_group

# if I set scale. = FALSE I get a downfacing horseshoe
df4.pca <- prcomp(clr_norm_df, center = TRUE,scale. = TRUE)
summary(df4.pca)

substr(rownames(clr_norm_df),1,3)

pdf("cm_PCA.pdf")
ggbiplot(df4.pca,labels=rownames(clr_norm_df),groups=substr(rownames(clr_norm_df),1,3),ellipse=TRUE,choices = (1:2))+
  theme_minimal()
dev.off()


#################################





