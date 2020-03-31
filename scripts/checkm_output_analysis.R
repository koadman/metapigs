# checkm_output of bins
# how many in total? 
# how many >90% and <5% contamination? 
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

######################################################################

# upload input 
concatenated_checkm_nearly <- read_delim("/Users/12705859/Desktop/metapigs_dry/checkm/checkm_all_nearly", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)

NROW(concatenated_checkm_nearly)
# some formatting 
concatenated_checkm_nearly$Completeness <- as.numeric(concatenated_checkm_nearly$Completeness)
concatenated_checkm_nearly$Contamination <- as.numeric(concatenated_checkm_nearly$Contamination)
# remove rows containing NAs as in the original file these rows where headers
concatenated_checkm_nearly <- na.omit(concatenated_checkm_nearly)

######################################################################

# -- division of nearly complete bins in quality categories ( 90&5 ; 95&5 ; 99&0.001 )

concatenated_checkm_nearly <- dplyr::filter(concatenated_checkm_nearly, !grepl("Completeness",Completeness))


newdata <- subset(concatenated_checkm_nearly, Completeness >= 90 & Contamination <= 5)
paste0((round(NROW(newdata)/NROW(concatenated_checkm_nearly),4)*100),"%")
nrow(newdata)

newdata1 <- subset(concatenated_checkm_nearly, Completeness >= 95 & Contamination <= 5)
paste0((round(NROW(newdata1)/NROW(concatenated_checkm_nearly),4)*100),"%")

newdata2 <- subset(concatenated_checkm_nearly, Completeness >= 99 & Contamination <= 0.001)
paste0((round(NROW(newdata2)/NROW(concatenated_checkm_nearly),4)*100),"%")

######################################################################

# -- proportion of >=90 and <=5 bins assigned to phylum/order/etc... : 

# how many of >90 <5 bins are classified at the each level? 
nn <- newdata %>% rename(taxa = `Taxonomy (contained)`)

archaea <- filter(nn, grepl('k__Archaea', taxa))
bacteria <- filter(nn, grepl('k__Bacteria', taxa))

kingdom <- filter(nn, grepl('k__', taxa))
NROW(kingdom)

phylum <- filter(nn, grepl('p__', taxa))
paste0((round(NROW(phylum)/NROW(nn),4)*100),"%")

class <- filter(nn, grepl('c__', taxa))
paste0((round(NROW(class)/NROW(nn),4)*100),"%")

order <- filter(nn, grepl('o__', taxa))
paste0((round(NROW(order)/NROW(nn),4)*100),"%")

family <- filter(nn, grepl('f__', taxa))
paste0((round(NROW(family)/NROW(nn),4)*100),"%")

genus <- filter(nn, grepl('g__', taxa))
paste0((round(NROW(genus)/NROW(nn),4)*100),"%")

species <- filter(nn, grepl('s__', taxa))
paste0((round(NROW(species)/NROW(nn),4)*100),"%")

######################################################################

# -- Proportions of each phylum member in the >90 <5 bins: 

phy <- filter(nn, grepl('k__Bacteria;p__', taxa))
Actinobacteria <- filter(nn, grepl('Actinobacteria', taxa))
Bacteroidetes <- filter(nn, grepl('Bacteroidetes', taxa))
Chlamydiae <- filter(nn, grepl('Chlamydiae', taxa))
Firmicutes <- filter(nn, grepl('Firmicutes', taxa))
Proteobacteria <- filter(nn, grepl('Proteobacteria', taxa))
Spirochaetes <- filter(nn, grepl('Spirochaetes', taxa))
Synergistetes <- filter(nn, grepl('Synergistetes', taxa))
Tenericutes <- filter(nn, grepl('Tenericutes', taxa))
Verrucomicrobia <- filter(nn, grepl('Verrucomicrobia', taxa))

print(paste0("Of the >90% complete bins with <5% contamination,"," ", NROW(bacteria)/NROW(nn)*100,"%", " ", "were assigned to Bacteria and ",
             NROW(archaea), " ", "were assigned to Archaea"))

# Of the >90% complete bins with <5% contamination, 98.5117111255693% were assigned to Bacteria and 183 were assigned to Archaea
print(paste0("number of bacterial bins resolved at the phylum level:"," ", NROW(phy)))

# number of bacterial bins resolved at the phylum level: 11491
print(paste0("Percentage belonging to: ", "Actinobacteria:"," ", NROW(Actinobacteria)/NROW(phy)*100, " ; ",
      "Bacteroidetes:"," ", NROW(Bacteroidetes)/NROW(phy)*100, " ; ",
      "Chlamydiae:"," ", NROW(Chlamydiae)/NROW(phy)*100, " ; ",
      "Firmicutes:"," ", NROW(Firmicutes)/NROW(phy)*100, " ; ",
      "Proteobacteria:"," ", NROW(Proteobacteria)/NROW(phy)*100, " ; ",
      "Spirochaetes:"," ", NROW(Spirochaetes)/NROW(phy)*100, " ; ",
      "Synergistetes:"," ", NROW(Synergistetes)/NROW(phy)*100, " ; ",
      "Tenericutes:"," ", NROW(Tenericutes)/NROW(phy)*100, " ; ",
      "Verrucomicrobia:"," ", NROW(Verrucomicrobia)/NROW(phy)*100," "))
# Percentage belonging to: Actinobacteria: 6.70959881646506 ; Bacteroidetes: 13.6019493516665 ; 
# Chlamydiae: 0.783221651727439 ; Firmicutes: 65.4860325472109 ; Proteobacteria: 2.21912801322774 ; 
# Spirochaetes: 1.70568270820642 ; Synergistetes: 0.269776346706118 ; Tenericutes: 9.21590810199286 ; 
# Verrucomicrobia: 0.00870246279697154

######################################################################

# to the purpose of merging checkm info to bins:
# first upload bins with counts (sample-dereplicated- output of 7.R)

no_reps_all <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/no_reps_all.csv", 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)
NROW(no_reps_all)
head(no_reps_all)
View(no_reps_all)
mean(no_reps_all$value)
min(no_reps_all$value)
max(no_reps_all$value)
summary(no_reps_all$value)



#no_reps_all <- no_reps_all %>%
#  filter(!value==0)
hist(no_reps_all$value,breaks=10000,xlim = c(0, 150000))
hist(no_reps_all$value,breaks=100000,xlim = c(0, 2000))
# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)

# rename cols to matching colnames between dataframes to merge  
colnames(newdata)[colnames(newdata)=="pigid"] <- "pig"
colnames(newdata)[colnames(newdata)=="Bin Id"] <- "bin"
colnames(newdata)[colnames(newdata)=="Taxonomy (contained)"] <- "taxa"

newdata <- newdata %>%
  select(pig,bin,taxa)

NROW(newdata)
# 12296
length(unique(paste0(no_reps_all$pig,no_reps_all$bin)))
# 50787

# percentage of nearly complete bins 
NROW(newdata)/length(unique(paste0(no_reps_all$pig,no_reps_all$bin)))*100
# 24.21092

head(no_reps_all)
head(newdata)

df <- right_join(no_reps_all, newdata, by=c("pig","bin"))
head(df)
NROW(df)


length(unique(paste0(df$pig,df$bin)))
# yes it's correct, we still have 12296 pig bins (same as 12296 from checkm)

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
                                    "Dscour",
                                    "ColiGuard", 
                                    "Neomycin",
                                    "NeoD",
                                    "NeoC"))

df
##################################

# treemap PHYLA (piglet samples, mothers exlcuded)

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_2) %>%
  filter(!cohort=="Mothers")
for_treemap <- na.omit(for_treemap)

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

pdf("Desktop/metapigs_dry/checkm/treemap_phyla.pdf")
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="4 to 10 weeks old piglets - Phyla distribution - from nearly complete MAGs", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()

# treemap CLASS

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_3) %>%
  filter(!cohort=="Mothers")
for_treemap <- na.omit(for_treemap)

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

pdf("Desktop/metapigs_dry/checkm/treemap_class.pdf")
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()

# treemap ORDER

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_4) %>%
  filter(!cohort=="Mothers")
for_treemap <- na.omit(for_treemap)

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

pdf("Desktop/metapigs_dry/checkm/treemap_order.pdf")
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()

##################################

# parallel coordinates - PHYLA

# keep only rows that contain phyla info, discard others; also exclude mothers
for_parallel_coo <- df %>%
  dplyr::select(cohort,pig,bin,date,taxa_2,value) %>%
  filter(!cohort=="Mothers")
for_parallel_coo <- na.omit(for_parallel_coo)

# remove "p__" before phylum 
for_parallel_coo[5] <- lapply(
  for_parallel_coo[5], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

numbers <- setDT(for_parallel_coo)[, .(Freq = .N), by = .(taxa_2,date)]
numbers

pdf("Desktop/metapigs_dry/checkm/parallel_coordinates_phyla.pdf")
ggplot(numbers, aes(x=date, y=Freq, group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection point",
       y = "MAGs frequency")
dev.off()

head(for_parallel_coo)
NROW(for_parallel_coo)

# general time change - unrooted
summs_for_parallel_coo <- for_parallel_coo %>% group_by(taxa_2,date) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 

gen_unrooted <- ggplot(summs_for_parallel_coo, aes(x=date, y=mean, group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "unrooted PD - mean")
gen_unrooted

summs_for_parallel_coo <- for_parallel_coo %>% group_by(taxa_2,date) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 

pdf("~/Desktop/metapigs_dry/checkm/rel_ab_phyla.pdf")
ggplot(summs_for_parallel_coo, aes(fill=taxa_2, y=mean, x=date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

summs_for_parallel_coo <- for_parallel_coo %>% group_by(taxa_2,date,cohort) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 

pdf("Desktop/metapigs_dry/checkm/rel_ab_phyla_cohorts.pdf")
ggplot(summs_for_parallel_coo, aes(fill=taxa_2, y=mean, x=date)) + 
  geom_bar(position="fill", stat="identity")+
  facet_wrap(~cohort)
dev.off()

# higher proportion of Tenericutes in Neomycin from t5 to t10 compared to other cohorts?
# Euryarchaeota "disapper" from t5 to t10. 
# Actinobacteria slowly increase in proportion with time
# igher bacteroidetes proportion between t3 and t4 in Neo than other cohorts? 
# Synergistetes prettu much disappear after t1 in Control, Dscour and ColiGuard,
# whereas they disappera later (after t2) in Neo and NeoD. 
# In NeoC in t4 still high proportion before gone. 
# Proteobacterai from t4 on higher proportion in Dscour, low in any other cohort. 



##################################

# parallel coordinates - ORDER

# keep only rows that contain phyla info, discard others; also exclude mothers
for_parallel_coo <- df %>%
  dplyr::select(cohort,pig,bin,date,taxa_4) %>%
  filter(!cohort=="Mothers")
for_parallel_coo <- na.omit(for_parallel_coo)

# remove "p__" before phylum 
for_parallel_coo[5] <- lapply(
  for_parallel_coo[5], 
  gsub, 
  pattern = "o__", 
  replacement = "", 
  fixed = TRUE)

numbers <- setDT(for_parallel_coo)[, .(Freq = .N), by = .(taxa_4,date)]
numbers

pdf("Desktop/metapigs_dry/checkm/parallel_coordinates_order.pdf")
ggplot(numbers, aes(x=date, y=Freq, group=taxa_4, color=taxa_4)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection point",
       y = "MAGs frequency")
dev.off()


##################################


# does breed have an influence of phyla composition ? 

# load details (breed, line, bday, mothers)
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"
details <- read_excel(paste0(basedir, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'pig'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'stig'
colnames(details)[colnames(details) == '...8'] <- 'breed'
details$pig <- gsub("G","",details$pig)
details$pig <- gsub("T","",details$pig)

details <- details %>%
  dplyr::select(pig,BIRTH_DAY,LINE,breed,stig,nurse)

# join df with details: 
df <- right_join(for_parallel_coo,details)


summs_for_parallel_coo <- df %>% 
  dplyr::select(breed,pig,bin,date,taxa_2,value) %>%
  group_by(taxa_2,date,breed) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 

pdf("~/Desktop/metapigs_dry/checkm/rel_ab_phyla_breed.pdf")
ggplot(summs_for_parallel_coo, aes(fill=taxa_2, y=mean, x=date)) + 
  geom_bar(position="fill", stat="identity")+
  facet_wrap(~breed)
dev.off()




# keep only rows that contain phyla info, discard others; also exclude mothers
for_parallel_coo <- df %>%
  dplyr::select(cohort,pig,bin,date,taxa_2) %>%
  filter(!cohort=="Mothers")
for_parallel_coo <- na.omit(for_parallel_coo)

# remove "p__" before phylum 
for_parallel_coo[5] <- lapply(
  for_parallel_coo[5], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

numbers <- setDT(for_parallel_coo)[, .(Freq = .N), by = .(cohort,taxa_2,date)]
numbers

pdf("~/Desktop/metapigs_dry/checkm/parallel_coordinates_phyla_breed.pdf")
ggplot(numbers, aes(x=date, y=Freq, group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection point",
       y = "MAGs frequency")+
  facet_grid(~cohort)
dev.off()





##################################

# normalized for cohort size 
df2 <- df %>%
  group_by(cohort) %>%
  mutate(value=value/sum(value))

# kingdom
kingdom <- ggplot(data=df2, aes(x=cohort, y=value, fill=taxa_1)) +
  geom_bar(stat="identity") + 
  ggtitle("Kindgom distribution among cohorts") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Kingdom")
kingdom

# phylum
phylum <- ggplot(data=df2, aes(x=cohort, y=value, fill=taxa_2)) +
  geom_bar(stat="identity") + 
  ggtitle("Phyla distribution among cohorts") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Phylum")
phylum

# class
class <- ggplot(data=df2, aes(x=cohort, y=value, fill=taxa_3)) +
  geom_bar(stat="identity") + 
  ggtitle("Class distribution among cohorts") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Class")
class

# order
order <- ggplot(data=df2, aes(x=cohort, y=value, fill=taxa_4)) +
  geom_bar(stat="identity") + 
  ggtitle("Order distribution among cohorts") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Order") +
  theme(legend.position="right") + guides(fill = guide_legend(ncol = 1))
order

pdf("/Users/12705859/Desktop/metapigs_dry/checkm/cohorts_kpco_distribution.pdf")
kingdom
phylum
class
order
dev.off()


##########################################################







# Barplots: phylum over time and cohorts (normalization by lib size)

ctrl <- df %>%
  filter(cohort == "Control") %>%
  group_by(pig) %>%
  mutate(value=value/sum(value)) %>%
  group_by(date) %>%
  mutate(value=value/sum(value))

Dscour <- df %>%
  filter(cohort == "Dscour") %>%
  group_by(pig) %>%
  mutate(value=value/sum(value)) %>%
  group_by(date) %>%
  mutate(value=value/sum(value))

ColiGuard <- df %>%
  filter(cohort == "ColiGuard") %>%
  group_by(pig) %>%
  mutate(value=value/sum(value)) %>%
  group_by(date) %>%
  mutate(value=value/sum(value))

neo <- df %>%
  filter(cohort == "Neomycin") %>%
  group_by(pig) %>%
  mutate(value=value/sum(value)) %>%
  group_by(date) %>%
  mutate(value=value/sum(value))

neoD <- df %>%
  filter(cohort == "NeoD") %>%
  group_by(pig) %>%
  mutate(value=value/sum(value)) %>%
  group_by(date) %>%
  mutate(value=value/sum(value))

neoC <- df %>%
  filter(cohort == "NeoC") %>%
  group_by(pig) %>%
  mutate(value=value/sum(value)) %>%
  group_by(date) %>%
  mutate(value=value/sum(value))

phylum_ctrl <- ggplot(data=ctrl, aes(x=date, y=value, fill=taxa_2)) +
  geom_bar(stat="identity") + 
  ggtitle("Phyla distribution over time - Control") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Phylum")
phylum_ctrl

phylum_Dscour <- ggplot(data=Dscour, aes(x=date, y=value, fill=taxa_2)) +
  geom_bar(stat="identity") + 
  ggtitle("Phyla distribution over time - D-scour") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Phylum")
phylum_Dscour

phylum_ColiGuard <- ggplot(data=ColiGuard, aes(x=date, y=value, fill=taxa_2)) +
  geom_bar(stat="identity") + 
  ggtitle("Phyla distribution over time - ColiGuard") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Phylum")
phylum_ColiGuard

phylum_neo <- ggplot(data=neo, aes(x=date, y=value, fill=taxa_2)) +
  geom_bar(stat="identity") + 
  ggtitle("Phyla distribution over time - Neomycin") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Phylum")
phylum_neo

phylum_neoD <- ggplot(data=neoD, aes(x=date, y=value, fill=taxa_2)) +
  geom_bar(stat="identity") + 
  ggtitle("Phyla distribution over time - Neomycin+D-scour") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Phylum")
phylum_neoD

phylum_neoC <- ggplot(data=neoC, aes(x=date, y=value, fill=taxa_2)) +
  geom_bar(stat="identity") + 
  ggtitle("Phyla distribution over time - Neomycin+ColiGuard") +
  ylab("ratio of reads mapping against nearly complete genomes") +
  labs(fill = "Phylum")
phylum_neoC


pdf("~/Desktop/metapigs_dry/checkm/phyla_distribution_per_cohort.pdf")
phylum_ctrl
phylum_Dscour
phylum_ColiGuard
phylum_neo
phylum_neoD
phylum_neoC
dev.off()

####################################

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


pdf("/Users/12705859/Desktop/metapigs_dry/checkm/phylum_neo.pdf")
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


pdf("/Users/12705859/Desktop/metapigs_dry/checkm/phylum_ctrl.pdf")
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

pdf("/Users/12705859/Desktop/metapigs_dry/checkm/phyla_time_balloon.pdf")
balloonplot(t(m), main = "Phyla distribution from (CheckM) nearly complete bins over time", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()






# careful cause we don't have the pigID to distinguish which bins to which sample
all_checkm_output <- read_delim("~/Desktop/metapigs_dry/checkm/all_checkm_output.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
all_checkm_output <- dplyr::filter(all_checkm_output, !grepl("Completeness",Completeness))
all_checkm_output$Completeness <- as.numeric(all_checkm_output$Completeness)
all_checkm_output$Contamination <- as.numeric(all_checkm_output$Contamination)
NROW(all_checkm_output)
clean <- all_checkm_output

clean <- clean %>%
  filter(Completeness >= 70)
NROW(clean)
clean <- clean %>%
  filter(Contamination <= 10)
NROW(clean)

clean$type <- ifelse(clean$Completeness >=90 & clean$Contamination <=5, "high", "low")


# removing data less than <80 Complete to plot neatly
clean <- clean %>%
  filter(Completeness >= 80)
View(clean)
pdf("/Users/12705859/Desktop/metapigs_dry/checkm/checkm_all_completeness_contamination.pdf")
gg <- ggplot(clean, aes(x=Completeness, y=Contamination)) + 
  geom_point(aes(color=factor(type)),size=0.1,shape=21)+
  ylab("Contamination (%)")+
  xlab("Completeness (%)")+
  xlim(70,100)+
  theme_minimal()+
  theme(axis.title.x = element_text(),
        axis.title.y = element_text())
gg
dev.off()


clean$`N50 (scaffolds)` <- as.numeric(clean$`N50 (scaffolds)`)
clean$`N50 (contigs)` <- as.numeric(clean$`N50 (contigs)`)
clean$`# contigs` <- as.numeric(clean$`# contigs`)




clean2 <- clean

# log10 N50 scaffolds
clean2$`N50 (contigs)` <- log10(as.numeric(clean2$`N50 (contigs)`))
`N50 (contigs)` <- clean2$`N50 (contigs)`
hist(`N50 (contigs)`)
#ggplot(data = clean2) + geom_histogram(aes(x = `N50 (contigs)`), binwidth = 10000)


# log10 contigs
clean2$`# contigs` <- log10(as.numeric(clean2$`# contigs`))
`# contigs` <- clean2$`# contigs`
hist(`# contigs`)
#nolog
clean2$`# contigs` <- as.numeric(clean2$`# contigs`)
`# contigs` <- clean2$`# contigs`
hist(`# contigs`)

ggplot(data = clean2) + geom_histogram(aes(x = `# contigs`), binwidth = 0.1)

ggplot(data = clean2) + 
  geom_histogram(aes(x = `# contigs`), binwidth = 0.01) +
  scale_x_log10()

ggplot(data = clean2) + 
  geom_histogram(aes(x = `N50 (contigs)`), binwidth = 0.01) +
  scale_x_log10()






clean2
View(clean2)

`N50 (log10)` <- clean2$`N50 (scaffolds)`
hist(`N50 (log10)`)


az <- clean %>%
  select(`N50 (scaffolds)`, `Bin Id`)
bz <- unique(az)
cz <- data.frame(table(az$`N50 (scaffolds)`))
head(cz)
View(cz)
