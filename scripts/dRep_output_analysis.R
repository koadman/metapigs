
# dRep_output_analysis.R
# analysis of dRep output 

library(readxl)
library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(splitstackshape)
library(pheatmap)
library(ggpubr)
library(robCompositions)
library(ggbiplot)



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

no_reps_all$primary_cluster <- paste0(no_reps_all$secondary_cluster)
no_reps_all <- cSplit(no_reps_all,"primary_cluster","_")
no_reps_all$primary_cluster_2 <- NULL
colnames(no_reps_all)[colnames(no_reps_all)=="primary_cluster_1"] <- "primary_cluster"

######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv("/Users/12705859/Desktop/metapigs_dry/gtdbtk/gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################


# create text file to contain dRep text output

sink(file = "dRep_numbers.txt", 
     append = FALSE, type = c("output"))
sink()


###########################


sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Number of dRep-clustered bins: ", NROW(C1) )
paste0("of which primary clusters: ", length(unique(C1$primary_cluster)) )
paste0("of which secondary clusters ", length(unique(C1$secondary_cluster)) )
sink()



########################################################################################################


# Extent of agreement between dRep and GTDBTK classification: 


df <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))


# Primary clusters: 


b <- df %>%
  dplyr::group_by(primary_cluster) %>%
  #dplyr::filter(n()>300) %>%               # this is optional; comment out to look at all
  dplyr::filter(!primary_cluster=="no_cluster") %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::select(node,domain,phylum,order,family,genus,species,primary_cluster) 

NROW(unique(b$primary_cluster))

# most frequent primary clusters 


b_species <- b %>%
  dplyr::group_by(primary_cluster,species) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2,species) %>%
  dplyr::rename(., species_agree = num2) 

b_genus <- b %>%
  dplyr::group_by(primary_cluster,genus) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2) %>%
  dplyr::rename(., genus_agree = num2) 

b_family <- b %>%
  dplyr::group_by(primary_cluster,family) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2,family) %>%
  dplyr::rename(., family_agree = num2) 

b_order <- b %>%
  dplyr::group_by(primary_cluster,order) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2) %>%
  dplyr::rename(., order_agree = num2) 


# saving these as variables before removing columns, to (optionally) use later as "fuller" taxa names
family_names <- b_family$family
b_family$family <- NULL
species_names <- b_species$species
b_species$species <- NULL

a <- merge(b_genus, b_species)
a <- merge(a,b_family)
a <- merge(a,b_order)
a

rownames(a) <- a[,1]
a[,1] <- NULL

summary_a <- as.data.frame.array(summary(a))


sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Extent of agreement between dRep classificantion and gtdbtk assignment of bins")
paste0("Primary clusters: ")
summary_a
sink()



######################################################################


# Secondary clusters: 


b <- df %>%
  dplyr::group_by(secondary_cluster) %>%
  #dplyr::filter(n()>300) %>%               # this is optional; comment out to look at all
  dplyr::filter(!secondary_cluster=="no_cluster") %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::select(node,domain,phylum,order,family,genus,species,secondary_cluster) 

NROW(unique(b$secondary_cluster))

# most frequent secondary clusters 


b_species <- b %>%
  group_by(secondary_cluster,species) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2,species) %>%
  dplyr::rename(., species_agree = num2) 

b_genus <- b %>%
  dplyr::group_by(secondary_cluster,genus) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2) %>%
  dplyr::rename(., genus_agree = num2) 

b_family <- b %>%
  dplyr::group_by(secondary_cluster,family) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2,family) %>%
  dplyr::rename(., family_agree = num2) 

b_order <- b %>%
  dplyr::group_by(secondary_cluster,order) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2) %>%
  dplyr::rename(., order_agree = num2) 


# saving these as variables before removing columns, to (optionally) use later as "fuller" taxa names
family_names <- b_family$family
b_family$family <- NULL
species_names <- b_species$species
b_species$species <- NULL

a <- merge(b_genus, b_species)
a <- merge(a,b_family)
a <- merge(a,b_order)
a

rownames(a) <- a[,1]
a[,1] <- NULL

summary_a <- as.data.frame.array(summary(a))


sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Extent of agreement between dRep classificantion and gtdbtk assignment of bins")
paste0("Secondary clusters: ")
summary_a
sink()




######################################################################
######################################################################

########################################################################################################



# Display amount of shared vs unique clusters 


# primary clusters, piglets
primary_piglets <- no_reps_all %>%
  filter(!cohort=="Mothers") %>%
  select(pig,primary_cluster) %>% 
  distinct() %>%
  group_by(primary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common primary clusters (95% ANI)", 
       subtitle="distribution among piglets",
       x = "piglets",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))

# secondary clusters, piglets
secondary_piglets <- no_reps_all %>%
  filter(!cohort=="Mothers") %>%
  select(pig,secondary_cluster) %>% 
  distinct() %>%
  group_by(secondary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common secondary clusters (99% ANI)", 
       subtitle="distribution among piglets",
       x = "piglets",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))




# primary clusters, mothers
primary_mothers <- no_reps_all %>%
  filter(cohort=="Mothers") %>%
  select(pig,primary_cluster) %>% 
  distinct() %>%
  group_by(primary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common primary clusters (95% ANI)", 
       subtitle="distribution among mothers",
       x = "mothers",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))

# secondary clusters, mothers
secondary_mothers <- no_reps_all %>%
  filter(cohort=="Mothers") %>%
  select(pig,secondary_cluster) %>% 
  distinct() %>%
  group_by(secondary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common secondary clusters (99% ANI)", 
       subtitle="distribution among mothers",
       x = "mothers",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))



tosave <- ggarrange(primary_piglets,
          secondary_piglets,
          primary_mothers,
          secondary_mothers,
          ncol=2,
          nrow=2,
          common.legend = TRUE,
          widths = c(1,1),
          heights = c(1,1))


pdf("dRep_common_vs_unique_clusters.pdf")
tosave
dev.off()


sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Percentage of shared primary clusters among piglets:")
NROW(which(primary_piglets$data$type=="common"))/NROW(primary_piglets$data)*100
paste0("Percentage of shared secondary clusters among piglets:")
NROW(which(secondary_piglets$data$type=="common"))/NROW(secondary_piglets$data)*100
paste0("#######")
paste0("Percentage of shared primary clusters among mothers:")
NROW(which(primary_mothers$data$type=="common"))/NROW(primary_mothers$data)*100
paste0("Percentage of shared secondary clusters among mothers:")
NROW(which(secondary_mothers$data$type=="common"))/NROW(secondary_mothers$data)*100
sink()



######################################################################################################
######################################################################################################

# Principal component analysis: clustering of bins based on phyla with time (labels per cohort)




# I hereby select taxa_2 only (corresponds to phylum) and remove any row where no phylum was resolved
df1 <- no_reps_all %>%
  select(pig,bin,date,value,secondary_cluster,cohort)

df1 <- as.data.frame(na.omit(df1))

unique(df1$secondary_cluster)


#################################
# STEP 1.

# normalization for library size 
df2 <- df1 %>%
  filter(!secondary_cluster=="no_cluster") %>%
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
  dplyr::group_by(pig,date,secondary_cluster) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

df3$sample<- paste0(df3$date,"_",df3$pig)
df3$pig <- NULL
df3$date <- NULL

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
  pivot_wider(names_from = secondary_cluster, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 
head(df4)
colnames(df4)
# # test:
# test3 <- test2 %>%
#   pivot_wider(names_from = secondary_cluster.y, values_from = indiv_sum, values_fill = list(indiv_sum = 0))
# head(test3)
# sum(test3[1,])

df4 <- as.data.frame(df4)
rowSums(df4[,-1])

#################################


# get a quick cohorts to pig table
cohorts <- df %>% dplyr::select(cohort,pig,date) %>% distinct()
cohorts$sample <- paste0(cohorts$date,"_",cohorts$pig)
cohorts <- as.data.frame(cohorts)


rownames(df4) <- df4$sample
df4$sample <- NULL
m <- as.matrix(df4)

df4.pca <- prcomp(m, center = FALSE,scale. = FALSE)
summary(df4.pca)

# to get samples info showing on PCA plot
this_mat_samples <- data.frame(sample=rownames(m))
this_mat_samples <- inner_join(this_mat_samples,cohorts)
NROW(this_mat_samples)


pdf("dRep_PCA.pdf")
ggbiplot(df4.pca,groups = this_mat_samples$date,ellipse=FALSE,var.axes = FALSE,choices = (1:2)) +
  theme_bw()
dev.off()


#################################




