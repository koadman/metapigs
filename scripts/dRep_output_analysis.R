# dRep.R
# analysis of dRep output 

library(readxl)
library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)




setwd("~/Desktop/metapigs_dry/dRep/")
basedir = "~/Desktop/metapigs_dry/"



# upload dRep output 
Cdb <- read_csv("Cdb.csv")

# upload cohorts info
cohorts <- read_xlsx(paste0(basedir,"cohorts.xlsx"))

C1 <- separate(data = Cdb, col = genome, into = c("pig", "bin"), sep = "_")
C1 <- C1[,c("pig","bin","primary_cluster","secondary_cluster")]
C1$primary_cluster <- as.character(C1$primary_cluster)
head(C1)


# upload bins with counts (from output of 7.R)
no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
                      na.strings=c("","NA"),
                      check.names = FALSE,
                      header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

##########################

# create text file to contain dRep text output

sink(file = "dRep_numbers.txt", 
     append = FALSE, type = c("output"))
sink()

###########################


# Bins distribution:
# how many bins each subject has:
az <- no_reps_all %>%
  dplyr::select(bin, pig, cohort) %>%
  distinct() 

az <- az %>%
  mutate(group = ifelse(cohort == "Mothers", "sows (n=42)","piglets (n=126)"))

az <- az %>%
  group_by(group,pig) %>%
  summarise(`number of metagenomes obtained`= n()) 
tail(az)


pdf("dRep_bins_distribution_subjects.pdf")
ggplot(data=az, mapping=aes(x=group, y=`number of metagenomes obtained`)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1)
dev.off()
  

# Get the numbers: 
az <- no_reps_all %>%
  dplyr::select(bin,pig)
bz <- unique(az)
cz <- data.frame(table(bz$pig))

sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("total of bins (excludes the pos and neg controls bins): ", sum(cz$Freq) )
paste0("Bins distribution across subjects: ")
summary(cz$Freq)
sink()


##########################

# Bins per timepoint: 

az <- no_reps_all %>%
  dplyr::select(bin, pig, date) %>%
  distinct() 


az <- az %>%
  group_by(date,pig) %>%
  summarise(`number of metagenomes obtained`= n()) 
tail(az)


pdf("dRep_bins_distribution_subjects.pdf")
ggplot(data=az, mapping=aes(x=date, y=`number of metagenomes obtained`)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1) +
  labs(title="Number of metagenomes obtained per timepoint", 
       x = "timepoint",
       y = "number of metagenomes obtained",
       subtitle=NULL)
dev.off()


##########################

# number of samples per time point: 

az <- no_reps_all %>%
  filter(!cohort=="Mothers") %>%
  dplyr::select(pig, date) %>%
  distinct() 

az <- az %>%
  group_by(date) %>%
  summarise(`number of samples`= n()) 
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

pdf("dRep_samples_per_subject.pdf")
ggplot(data=az, mapping=aes(x=date, y=`number of samples`, color=sampling)) + 
  geom_point() +
  geom_line(aes(group = sampling), linetype = 2) +
  theme_bw() + 
  labs(title="Number of samples per timepoint", 
                  x = "number of samples",
                  y = "timepoints",
                  subtitle=NULL)
dev.off()








# continue cleaning from here 




sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Timepoints per subject: ")
summary(cz2$Freq)
sink()

##########################

# Correlation between timepoints available and bins obtained per subject: 
# How many timepoints scale to how many bins? 

colnames(cz2)[colnames(cz2)=="Freq"] <- "Freq2"

bins_timepoints <- cbind(cz,cz2)
bins_timepoints$Var1 <- NULL

pdf("/Users/12705859/Desktop/metapigs_dry/dRep_timepoints_to_bins.pdf")
ggplot(data = bins_timepoints, aes(x = Freq, y = Freq2)) + 
  geom_point(color='blue') +
  geom_smooth(method = "lm", se = FALSE)+
  labs(title="Bins obtained versus number of timepoints available from each subject", 
       x = "number of bins obtained per subject",
       y = "number of timepoints available per subject",
       subtitle=NULL)
dev.off()


##########################

# Bins per cohort: 
az2 <- no_reps_all %>%
  select(cohort,bin,pig)
bz2 <- unique(az2)
cz2 <- data.frame(table(bz2$cohort))
sum(cz2$Freq)
pdf("dRep_bins_per_cohort.pdf")
g <- ggplot(cz2, aes(Var1, Freq))
g + geom_bar(stat="identity", width = 0.3, fill="tomato2") + 
  labs(title="Bins distribution across subjects", 
       x = "cohorts",
       y = "Frequency",
       subtitle=NULL) +
  theme(axis.title.x=element_text(),
        axis.text.x=element_text())
dev.off()

##########################


sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Number of dRep-clustered bins: ", NROW(C1) )
paste0("of which primary clusters: ", length(unique(C1$primary_cluster)) )
paste0("of which secondary clusters ", length(unique(C1$secondary_cluster)) )
sink()

########################################################################################################



# left join because there is clustered bins from dRep (Cdb) we 
# can t merge to no_reps_all as we excluded bins from the pos controls 
# in script 1_5_HPC.R


NROW(no_reps_all)
NROW(C1)
head(no_reps_all)
head(C1)

new_dataset <- no_reps_all %>% left_join(C1, by=c("pig","bin"))
head(new_dataset)
NROW(new_dataset)


# these are the frequencies of most frequent (>300) 
# unique primary clusters 
a <- new_dataset %>%
  group_by(primary_cluster) %>%
  filter(!is.na(primary_cluster)) %>%
  filter(n()>300)
a

counts <- setDT(new_dataset)[, .(Freq = .N), by = .(secondary_cluster.x,cohort)]
View(counts)




length(unique(a$primary_cluster))
pdf("/Users/12705859/Desktop/metapigs_dry/most_frequent_primary_clusters_distribution_cohorts.pdf")
g <- ggplot(a, aes(primary_cluster))
g + geom_bar(aes(fill=date), width = 0.5) + 
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text()) +
  labs(title="Distribution among cohorts of most common primary clusters (95% ANI)", 
       subtitle="Primary clusters appearing >300 times in the dataset",
       x="most common (n=190) clusters (95% ANI)",
       y="Frequency")
dev.off()


# these are the frequencies of most frequent (>300) 
# unique secondary clusters
a <- new_dataset %>%
  group_by(secondary_cluster) %>%
  filter(!is.na(secondary_cluster)) %>%
  filter(n()>300)
length(unique(a$secondary_cluster))
pdf("/Users/12705859/Desktop/metapigs_dry/most_frequent_secondary_clusters_distribution_cohorts.pdf")
g <- ggplot(a, aes(secondary_cluster))
g + geom_bar(aes(fill=cohort), width = 0.5) + 
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text()) +
  labs(title="Distribution among cohorts of most common secondary clusters (99% ANI)", 
       subtitle="Secondary clusters appearing >300 times in the dataset",
       x="most common (n=141) clusters (99% ANI)",
       y="Frequency")
dev.off()
# what species do these common secondary clusters correspond to ? 
# testing question in other script: dRep_to_taxa.R


###################################################

# Display amount of shared vs unique clusters 

# separate mothers from piglets 
a <- merge.data.frame(C1, cohorts, by.x="pig", by.y = "Animal ID")

mothers <- a %>% filter(
  `Cohort Name` == "Mothers")
piglets <- a %>% filter(
  `Cohort Name` == "Control" |
    `Cohort Name` == "Neomycin" |
    `Cohort Name` == "ColiGuard" |
    `Cohort Name` == "D-scour" |
    `Cohort Name` == "Neomycin+D-scour" |
    `Cohort Name` == "Neomycin+ColiGuard"
)


# in piglets
a <- piglets %>%
  select(pig,bin,primary_cluster,secondary_cluster)
# primary clusters: 
a <- a %>% 
  group_by(primary_cluster) %>% 
  mutate(type = ifelse(n() > 1, "common","unique"))
pdf("/Users/12705859/Desktop/metapigs_dry/dRep/shared_vs_common_primary_piglets.pdf")
g <- ggplot(a, aes(pig))
g + geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common primary clusters (95% ANI)", 
       subtitle="distribution among piglets",
       x = "piglets",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text())
dev.off()
length(which(a$type=="unique"))
NROW(a)
length(which(a$type=="unique"))/NROW(a)*100
# 291 of 21400 ( 1.36% ) unique primary clusters
# secondary clusters: 
a <- a %>% 
  group_by(secondary_cluster) %>% 
  mutate(type = ifelse(n() > 1, "common","unique"))
pdf("/Users/12705859/Desktop/metapigs_dry/dRep/shared_vs_common_secondary_piglets.pdf")
g <- ggplot(a, aes(pig))
g + geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common secondary clusters (99% ANI)", 
       subtitle="distribution among piglets",
       x = "piglets",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text())
dev.off()
length(which(a$type=="unique"))
NROW(a)
length(which(a$type=="unique"))/NROW(a)*100
# 2951 of 21400 ( 13.79% ) unique secondary clusters

# in mothers
a <- mothers %>%
  select(pig,bin,primary_cluster,secondary_cluster)
# primary clusters: 
a <- a %>% 
  group_by(primary_cluster) %>% 
  mutate(type = ifelse(n() > 1, "common","unique"))
pdf("/Users/12705859/Desktop/metapigs_dry/dRep/shared_vs_common_primary_mothers.pdf")
g <- ggplot(a, aes(pig))
g + geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common primary clusters (95% ANI)", 
       subtitle="distribution among mothers",
       x = "mothers",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text())
dev.off()
length(which(a$type=="unique"))
NROW(a)
length(which(a$type=="unique"))/NROW(a)*100
# 76 of 839 ( 9.06%% ) unique primary clusters
# secondary clusters: 
a <- a %>% 
  group_by(secondary_cluster) %>% 
  mutate(type = ifelse(n() > 1, "common","unique"))
pdf("/Users/12705859/Desktop/metapigs_dry/dRep/shared_vs_common_secondary_mothers.pdf")
g <- ggplot(a, aes(pig))
g + geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common secondary clusters (99% ANI)", 
       subtitle="distribution among mothers",
       x = "mothers",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text())
dev.off()
length(which(a$type=="unique"))
NROW(a)
length(which(a$type=="unique"))/NROW(a)*100
# 246 of 839 ( 29.32% ) unique primary clusters













######################################################################################################
######################################################################################################

# Principal component analysis: clustering of bins based on phyla with time (labels per cohort)




# I hereby select taxa_2 only (corresponds to phylum) and remove any row where no phylum was resolved
df1 <- new_dataset %>%
  select(pig,bin,date,value,secondary_cluster.y,cohort)

df1 <- as.data.frame(na.omit(df1))

unique(df1$secondary_cluster.y)


#################################
# STEP 1.

# normalization for library size 
df2 <- df1 %>%
  dplyr::group_by(pig,date) %>%
  dplyr::mutate(norm_value=value/sum(value)) 
NROW(df1)
head(df1)

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
  dplyr::group_by(pig,date,secondary_cluster.y) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

# # test:
# test2 <- test %>%
#   group_by(secondary_cluster.y) %>%
#   dplyr::summarise(indiv_sum = sum(norm_value))
# head(test2)
# sum(test2$indiv_sum)

#################################
# STEP 3.

# long to wide format
df4 <- df3 %>%
  pivot_wider(names_from = secondary_cluster.y, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 
head(df4)

# # test:
# test3 <- test2 %>%
#   pivot_wider(names_from = secondary_cluster.y, values_from = indiv_sum, values_fill = list(indiv_sum = 0))
# head(test3)
# sum(test3[1,])


#################################

# STEP 4. PCA, dots are cohort_date 

# get a quick cohorts to pig table 
cohorts <- new_dataset %>% dplyr::select(cohort,pig) %>% distinct()

# join the cohort info
df5 <- inner_join(df4,cohorts) %>%
  dplyr::mutate(coho_date_group=paste0(date,"_",cohort)) 
df5
colnames(df5)
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
ggbiplot(df4.pca,labels=rownames(clr_norm_df),groups=substr(rownames(clr_norm_df),1,3),ellipse=TRUE,choices = (1:2))
dev.off()


#################################




