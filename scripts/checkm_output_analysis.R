

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

basedir = "/Users/12705859/Desktop/metapigs_dry/"
setwd("/Users/12705859/Desktop/metapigs_dry/checkm/")

# input files: 
# all_checkm_output.tsv
# checkm/checkm_all_nearly
# no_reps_all.csv
# pigTrial_GrowthWtsGE.hlsx.xlsx


# output files: 
# no_reps_all.csv
# 

######################################################################

# template to collect and store checkM info 

sink("checkm_numbers.txt")
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

pdf("checkm_Compl_vs_Contam.pdf")
ggplot(df_80_10, aes(x=Completeness, y=Contamination)) + 
  geom_point(aes(color=factor(type)),size=0.1,shape=21)+
  ylab("Contamination (%)")+
  xlab("Completeness (%)")+
  xlim(70,100)+
  theme_minimal()+
  theme(axis.title.x = element_text(),
        axis.title.y = element_text())
dev.off()

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

sink(file = "checkm_numbers.txt", 
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


pdf("checkm_contigs_distribution.pdf")
# Distribution of contig number across bins 
hist(df_90_5$`# contigs`, 
     main = "Distribution of # of contigs across nearly complete bins (>=90% <=5%)",
     breaks=100,
     xlab = "contigs")
hist(log10(df_90_5$`# contigs`), 
     main = "Distribution of # of contigs across nearly complete bins (>=90% <=5%) - log scale",
     breaks=100,
     xlab = "log10(contigs)")
hist(df_90_5$`N50 (scaffolds)`, 
     main = "Distribution of N50 (scaffolds) across nearly complete bins (>=90% <=5%)",
     breaks=100,
     xlab = "log10(N50 scaffolds)",
     xlim=c(0,3e+05))
hist(log10(df_90_5$`N50 (scaffolds)`), 
     main = "Distribution of N50 (scaffolds) across nearly complete bins (>=90% <=5%) - log scale",
     breaks=100,
     xlab = "log10(N50 scaffolds)",
     xlim=c(3.5,6))
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

sink(file = "checkm_numbers.txt", 
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

# -- Proportions of each phylum member in the >90 <5 bins: 

phy <- filter(nn, grepl('p__', taxa))
Actinobacteria <- filter(nn, grepl('Actinobacteria', taxa))
Bacteroidetes <- filter(nn, grepl('Bacteroidetes', taxa))
Chlamydiae <- filter(nn, grepl('Chlamydiae', taxa))
Euryarchaeota <- filter(nn, grepl('Euryarchaeota', taxa))
Firmicutes <- filter(nn, grepl('Firmicutes', taxa))
Proteobacteria <- filter(nn, grepl('Proteobacteria', taxa))
Spirochaetes <- filter(nn, grepl('Spirochaetes', taxa))
Synergistetes <- filter(nn, grepl('Synergistetes', taxa))
Tenericutes <- filter(nn, grepl('Tenericutes', taxa))
Verrucomicrobia <- filter(nn, grepl('Verrucomicrobia', taxa))



sink(file = "checkm_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("######################################################################")
paste0(" Phyla distribution : ")
paste0("Phylum Actinobacteria:   ",
       NROW(Actinobacteria),"  ",round(NROW(Actinobacteria)/NROW(phy),4)*100,"%")
paste0("Phylum Bacteroidetes:   ",
       NROW(Bacteroidetes),"  ",round(NROW(Bacteroidetes)/NROW(phy),4)*100,"%")
paste0("Phylum Chlamydiae:   ",
       NROW(Chlamydiae),"  ",round(NROW(Chlamydiae)/NROW(phy),4)*100,"%")
paste0("Phylum Euryarchaeota:   ",
       NROW(Euryarchaeota),"  ",round(NROW(Euryarchaeota)/NROW(phy),4)*100,"%")
paste0("Phylum Firmicutes:   ",
       NROW(Firmicutes),"  ",round(NROW(Firmicutes)/NROW(phy),4)*100,"%")
paste0("Phylum Proteobacteria:   ",
       NROW(Proteobacteria),"  ",round(NROW(Proteobacteria)/NROW(phy),4)*100,"%")
paste0("Phylum Spirochaetes:   ",
       NROW(Spirochaetes),"  ",round(NROW(Spirochaetes)/NROW(phy),4)*100,"%")
paste0("Phylum Synergistetes:   ",
       NROW(Synergistetes),"  ",round(NROW(Synergistetes)/NROW(phy),4)*100,"%")
paste0("Phylum Tenericutes:   ",
       NROW(Tenericutes),"  ",round(NROW(Tenericutes)/NROW(phy),4)*100,"%")
paste0("Phylum Verrucomicrobia:   ",
       NROW(Verrucomicrobia),"  ",round(NROW(Verrucomicrobia)/NROW(phy),4)*100,"%")
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

pdf("treemap_phyla.pdf")
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

pdf("treemap_class.pdf")
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

pdf("treemap_order.pdf")
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

# general time change - unrooted
summs_for_parallel_coo <- for_parallel_coo %>% group_by(taxa_2,date) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 

pdf("parallel_coordinates_phyla.pdf")
ggplot(summs_for_parallel_coo, aes(x=date, y=log10(mean), group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "log10 of the mean",
       color = "Phylum")  +
  theme(legend.title = element_text()) 
dev.off()


pdf("rel_ab_phyla.pdf")
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

pdf("rel_ab_phyla_cohorts.pdf")
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

library(robCompositions)


summs_for_parallel_coo <- for_parallel_coo %>% group_by(taxa_2,date) %>% 
  dplyr::summarise(value = sum(value))


test <- summs_for_parallel_coo %>%
  dplyr::select(taxa_2,date,value) %>%
  pivot_wider(names_from = taxa_2, values_from = value, values_fill = list(value = 0))

test[,11] <- NULL
labels <- test[,1]

test <- cenLR(test[,-1])
test <- test$x.clr

test <- as.data.frame(test)
test <- cbind(test,labels)
rownames(test) <- test[,10]
test[,10] <- NULL
test <- as.table(as.matrix(test))

balloonplot(t(test), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)



# Lib size normalization
# 
# Balloon plots per cohort

# lib size normalization
df1 <- df %>%
  group_by(pig,date) %>%
  mutate(value=value/sum(value)) %>% 
  select(pig, date, cohort,taxa_2, value)

# omit bins that don't contain phyla info
df1 <- na.omit(df1)

# aggregate: mean of pig-bin-date values from same cohort
x <- df1 %>%
  group_by(cohort, date, taxa_2) %>%
  summarise_each(funs(mean=mean(., na.rm = TRUE)))
x$pig_mean <- NULL
# now we have one value for each phylum-cohort-date
length(unique(paste0(x$cohort,x$date)))
# it's 6 cohorts times 10 timepoints, correct. 

# each cohort+date is a group : 
x$group <- paste0(x$cohort,"_",x$date) 
x$cohort <- NULL
x$date <- NULL
# remove the p__ suffix
x$taxa_2 <- gsub("p__","",x$taxa_2)
x

# normalize by group 
x <- x %>% 
  group_by(group) %>% mutate(value_mean = value_mean/sum(value_mean)*100)

View(x)
x2 <- x %>%
  pivot_wider(names_from = taxa_2, values_from = value_mean)

x2[is.na(x2)] <- 0

# selection of cohorts to plot individual balloon plots (1 per cohort)
# cohorts selection
x3_1 <- x2[grepl("Neomycin", x2[["group"]]), ]
rownames(x3_1) <- x3_1$group
x3_1$group <- NULL
dt_1 <- as.table(as.matrix(x3_1))
# cohorts selection
x3_2 <- x2[grepl("NeoD", x2[["group"]]), ]
rownames(x3_2) <- x3_2$group
x3_2$group <- NULL
dt_2 <- as.table(as.matrix(x3_2))
# cohorts selection
x3_3 <- x2[grepl("NeoC", x2[["group"]]), ]
rownames(x3_3) <- x3_3$group
x3_3$group <- NULL
x3_3[is.na(x3_3)] <- 0
dt_3 <- as.table(as.matrix(x3_3))

# cohorts selection
x3_4 <- x2[grepl("Control", x2[["group"]]), ]
rownames(x3_4) <- x3_4$group
x3_4$group <- NULL
dt_4 <- as.table(as.matrix(x3_4))
# cohorts selection
x3_5 <- x2[grepl("Dscour", x2[["group"]]), ]
rownames(x3_5) <- x3_5$group
x3_5$group <- NULL
dt_5 <- as.table(as.matrix(x3_5))
# cohorts selection
x3_6 <- x2[grepl("ColiGuard", x2[["group"]]), ]
rownames(x3_6) <- x3_6$group
x3_6$group <- NULL
x3_6[is.na(x3_6)] <- 0
dt_6 <- as.table(as.matrix(x3_6))


pdf("phylum_neo.pdf")
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
balloonplot(t(dt_1), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(dt_2), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(dt_3), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)

dev.off()


pdf("phylum_ctrl.pdf")
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
balloonplot(t(dt_4), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(dt_5), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(dt_6), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()





####################################

#balloon plot all cohorts - time

NROW(unique(paste0(df$pig,df$bin)))
head(df)

# normalization for library size 
df1 <- df %>%
  group_by(pig,date) %>%
  mutate(norm_value=value/sum(value)) %>% 
  select(pig, date, taxa_2, value, norm_value) 

# sum all the norm values that fall within same pig,date,taxa_2
df2 <- df1 %>%
  group_by(pig,date,taxa_2) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
  
# take the mean of each taxa by date 
df3 <- df2 %>%
  group_by(taxa_2,date) %>%
  dplyr::summarize(mean = mean(indiv_sum, na.rm=TRUE)) %>%
  mutate(perc=mean*100) %>%
  select(taxa_2,date,perc)

# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(x)


df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc)

df4 <- as.data.frame(df4)
rownames(df4) <- df4$date
df4$date <- NULL

m <- as.table(as.matrix(df4))

pdf("phyla_time_balloon.pdf")
balloonplot(t(m), main = "Phyla distribution from (CheckM) nearly complete bins over time", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()







