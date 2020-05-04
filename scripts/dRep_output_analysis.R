
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



######################################################################


# create text file to contain dRep text output

sink(file = "dRep_numbers.txt", 
     append = FALSE, type = c("output"))
sink()


###########################


sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("dRep-clustered bins: ", 
       round(NROW(C1)/NROW(gtdbtk_bins)*100,2),
       "%",
       " (n=",NROW(C1),")" )
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
  dplyr::select(node,domain,phylum,class,order,family,genus,species,primary_cluster) 

b_species <- b %>%
  dplyr::group_by(primary_cluster,species) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2,species) %>%
  dplyr::rename(., species_agree = num2) 
b_species <- as.data.frame.array(summary(b_species))
b_species <- as.data.frame(b_species[,2])

b_genus <- b %>%
  dplyr::group_by(primary_cluster,genus) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2) %>%
  dplyr::rename(., genus_agree = num2) 
b_genus <- as.data.frame.array(summary(b_genus))
b_genus <- as.data.frame(b_genus[,2])

b_family <- b %>%
  dplyr::group_by(primary_cluster,family) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2,family) %>%
  dplyr::rename(., family_agree = num2) 
b_family <- as.data.frame.array(summary(b_family))
b_family <- as.data.frame(b_family[,2])

b_order <- b %>%
  dplyr::group_by(primary_cluster,order) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2) %>%
  dplyr::rename(., order_agree = num2) 
b_order <- as.data.frame.array(summary(b_order))
b_order <- as.data.frame(b_order[,2])

b_class <- b %>%
  dplyr::group_by(primary_cluster,class) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2) %>%
  dplyr::rename(., class_agree = num2) 
b_class <- as.data.frame.array(summary(b_class))
b_class <- as.data.frame(b_class[,2])

b_phylum <- b %>%
  dplyr::group_by(primary_cluster,phylum) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,num2) %>%
  dplyr::rename(., phylum_agree = num2) 
b_phylum <- as.data.frame.array(summary(b_phylum))
b_phylum <- as.data.frame(b_phylum[,2])

prim_clu_agree <- cbind(b_phylum,
      b_class,
      b_order,
      b_family,
      b_genus,
      b_species)
colnames(prim_clu_agree) <- c("phylum","class",
                             "order","family",
                             "genus","species")

sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Extent of agreement between dRep classification and gtdbtk assignment of bins")
paste0("Primary clusters: ")
prim_clu_agree
sink()


######################################################################


# Secondary clusters: 


b <- df %>%
  dplyr::group_by(secondary_cluster) %>%
  #dplyr::filter(n()>300) %>%               # this is optional; comment out to look at all
  dplyr::filter(!secondary_cluster=="no_cluster") %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::select(node,domain,phylum,class,order,family,genus,species,secondary_cluster) 


b_species <- b %>%
  dplyr::group_by(secondary_cluster,species) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2,species) %>%
  dplyr::rename(., species_agree = num2) 
b_species <- as.data.frame.array(summary(b_species))
b_species <- as.data.frame(b_species[,2])

b_genus <- b %>%
  dplyr::group_by(secondary_cluster,genus) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2) %>%
  dplyr::rename(., genus_agree = num2) 
b_genus <- as.data.frame.array(summary(b_genus))
b_genus <- as.data.frame(b_genus[,2])

b_family <- b %>%
  dplyr::group_by(secondary_cluster,family) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2,family) %>%
  dplyr::rename(., family_agree = num2) 
b_family <- as.data.frame.array(summary(b_family))
b_family <- as.data.frame(b_family[,2])

b_order <- b %>%
  dplyr::group_by(secondary_cluster,order) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2) %>%
  dplyr::rename(., order_agree = num2) 
b_order <- as.data.frame.array(summary(b_order))
b_order <- as.data.frame(b_order[,2])

b_class <- b %>%
  dplyr::group_by(secondary_cluster,class) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2) %>%
  dplyr::rename(., class_agree = num2) 
b_class <- as.data.frame.array(summary(b_class))
b_class <- as.data.frame(b_class[,2])

b_phylum <- b %>%
  dplyr::group_by(secondary_cluster,phylum) %>%
  dplyr::summarise(num= n())  %>%
  dplyr::mutate(num2= num/sum(num))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,num2) %>%
  dplyr::rename(., phylum_agree = num2) 
b_phylum <- as.data.frame.array(summary(b_phylum))
b_phylum <- as.data.frame(b_phylum[,2])
class(sec_clu_agree)
sec_clu_agree <- cbind(b_phylum,
                        b_class,
                        b_order,
                        b_family,
                        b_genus,
                        b_species)
colnames(sec_clu_agree) <- c("phylum","class",
                             "order","family",
                             "genus","species")

sink(file = "dRep_numbers.txt", 
     append = TRUE, type = c("output"))
paste0("Extent of agreement between dRep classification and gtdbtk assignment of bins")
paste0("Secondary clusters: ")
sec_clu_agree
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

dRep_PC12 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (1:2)) +
  theme_bw() +
  xlim(c(-2,1)) +
  guides(color = guide_legend(nrow = 1))
dRep_PC34 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (3:4)) +
  theme_bw() +
  guides(color = guide_legend(nrow = 1))


dRep_PCA <- ggarrange(dRep_PC12,dRep_PC34,ncol=2,
                      common.legend=TRUE)

pdf("dRep_PCA.pdf")
annotate_figure(dRep_PCA,
                top = text_grob("PCA from piglets' cluster-assigned MAGs (dRep) (n=22403)",
                                size = 13))
dev.off()

