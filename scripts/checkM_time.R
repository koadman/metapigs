
library(pheatmap)
library(ggbiplot)
library(robCompositions)
library(splitstackshape)


setwd("/Users/12705859/Desktop/metapigs_dry/checkm/")


df <- right_join(no_reps_all, newdata, by=c("pig","bin"))
head(df)
NROW(df)

# split taxa column into several (kingdom, phylum, etc ...) 
df <- cSplit(df, "taxa", sep=";")
unique(df$cohort)

# keep only piglets cohorts
df <- df %>%
  filter(cohort=="Control"| 
         cohort=="Dscour"|
         cohort=="ColiGuard"| 
         cohort=="Neomycin"|
         cohort=="NeoD"|
         cohort=="NeoC")

# reorder dates 
df$date  = factor(df$date, levels=c("t0",
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

unique(df$cohort)
head(df)

# I hereby select taxa_2 only (corresponds to phylum) and remove any row where no phylum was resolved 
df <- df %>%
  select(pig,bin,date,value,taxa_2,cohort) %>%
  na.omit(.)

# remove "p__" before phylum 
df[5] <- lapply(
  df[5], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

unique(df$taxa_2)

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


# test <- df3 %>%
#   filter(date=="t6")
# sum(test$perc)


# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(df3)

df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

# prepare as matrix to plot
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



####################################

#balloon plot time - CONTROL
unique(df$cohort)
df0 <- df %>%
  filter(cohort=="Control") %>%
  select(pig,bin,date,value,taxa_2)
head(df0)
NROW(df0)

# normalization for library size 
df1 <- df0 %>%
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


# test <- df3 %>%
#   filter(date=="t6")
# sum(test$perc)


# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(df3)

df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

# prepare as matrix to plot
df4 <- as.data.frame(df4)
rownames(df4) <- df4$date
df4$date <- NULL
m_control <- as.table(as.matrix(df4))


####################################

#balloon plot time - Dscour
unique(df$cohort)
df0 <- df %>%
  filter(cohort=="Dscour") %>%
  select(pig,bin,date,value,taxa_2)
head(df0)
NROW(df0)

# normalization for library size 
df1 <- df0 %>%
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


# test <- df3 %>%
#   filter(date=="t6")
# sum(test$perc)


# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(df3)

df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

# prepare as matrix to plot
df4 <- as.data.frame(df4)
rownames(df4) <- df4$date
df4$date <- NULL
m_Dscour <- as.table(as.matrix(df4))

####################################

#balloon plot time - ColiGuard
unique(df$cohort)
df0 <- df %>%
  filter(cohort=="ColiGuard") %>%
  select(pig,bin,date,value,taxa_2)
head(df0)
NROW(df0)

# normalization for library size 
df1 <- df0 %>%
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


# test <- df3 %>%
#   filter(date=="t6")
# sum(test$perc)


# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(df3)

df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

# prepare as matrix to plot
df4 <- as.data.frame(df4)
rownames(df4) <- df4$date
df4$date <- NULL
m_ColiGuard <- as.table(as.matrix(df4))


####################################

#balloon plot time - Neomycin
unique(df$cohort)
df0 <- df %>%
  filter(cohort=="Neomycin") %>%
  select(pig,bin,date,value,taxa_2)
head(df0)
NROW(df0)

# normalization for library size 
df1 <- df0 %>%
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


# test <- df3 %>%
#   filter(date=="t6")
# sum(test$perc)


# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(df3)

df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

# prepare as matrix to plot
df4 <- as.data.frame(df4)
rownames(df4) <- df4$date
df4$date <- NULL
m_Neomycin <- as.table(as.matrix(df4))

####################################

#balloon plot time - NeoD
unique(df$cohort)
df0 <- df %>%
  filter(cohort=="NeoD") %>%
  select(pig,bin,date,value,taxa_2)
head(df0)
NROW(df0)

# normalization for library size 
df1 <- df0 %>%
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


# test <- df3 %>%
#   filter(date=="t6")
# sum(test$perc)


# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(df3)

df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

# prepare as matrix to plot
df4 <- as.data.frame(df4)
rownames(df4) <- df4$date
df4$date <- NULL
m_NeoD <- as.table(as.matrix(df4))


####################################

#balloon plot time - NeoC
unique(df$cohort)
df0 <- df %>%
  filter(cohort=="NeoC") %>%
  select(pig,bin,date,value,taxa_2)
head(df0)
NROW(df0)

# normalization for library size 
df1 <- df0 %>%
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


# test <- df3 %>%
#   filter(date=="t6")
# sum(test$perc)


# now we have one value for each phylum-cohort-date (initial df same as last df)
length(unique(paste0(df$taxa_2,df$date))) == NROW(df3)

df4 <- df3 %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

# prepare as matrix to plot
df4 <- as.data.frame(df4)
rownames(df4) <- df4$date
df4$date <- NULL
m_NeoC <- as.table(as.matrix(df4))


pdf("/Users/12705859/Desktop/metapigs_dry/checkm/phyla_time_balloon_bycohort.pdf")
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
balloonplot(t(m_control), main = "Phyla distribution from (CheckM) nearly complete bins over time - Control", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(m_Dscour), main = "Phyla distribution from (CheckM) nearly complete bins over time - D-scour", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(m_ColiGuard), main = "Phyla distribution from (CheckM) nearly complete bins over time - ColiGuard", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
balloonplot(t(m_Neomycin), main = "Phyla distribution from (CheckM) nearly complete bins over time - Neomycin", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(m_NeoD), main = "Phyla distribution from (CheckM) nearly complete bins over time - Neomycin+D-scour", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
balloonplot(t(m_NeoC), main = "Phyla distribution from (CheckM) nearly complete bins over time - Neomycin+ColiGuard", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()


pheatmap(t(m_control), display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
pheatmap(t(m_Dscour), display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
pheatmap(t(m_ColiGuard), display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
pheatmap(t(m_Neomycin), display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
pheatmap(t(m_NeoD), display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
pheatmap(t(m_NeoC), display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)


##################################################################

# Principal component analysis: clustering of bins based on phyla with time (labels per cohort)

#################################
# STEP 1.

# normalization for library size 
df1 <- df %>%
  dplyr::group_by(pig,date) %>%
  dplyr::mutate(norm_value=value/sum(value)) 
NROW(df1)
head(df1)

# # test:
# test <- df %>%
#   filter(pig=="14159") %>%
#   filter(date=="t0") %>%
#   dplyr::mutate(norm_value=value/sum(value))
# head(test)
# sum(test$norm_value)

#################################
# STEP 2.

# sum all the norm values that fall within same pig,date,taxa_2
df2 <- df1 %>%
  dplyr::group_by(pig,date,taxa_2) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df2)

# # test: 
# test2 <- test %>%
#   group_by(taxa_2) %>%
#   dplyr::summarise(indiv_sum = sum(norm_value))
# head(test2)
# sum(test2$indiv_sum)

#################################
# STEP 3.

# long to wide format
df3 <- df2 %>%
  pivot_wider(names_from = taxa_2, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 
head(df3)

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
df5 <- inner_join(df3,cohorts) %>%
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

pdf("checkM_phyla_time_PCA.pdf")
ggbiplot(df4.pca,labels=rownames(clr_norm_df),groups=substr(rownames(clr_norm_df),1,3),ellipse=TRUE,choices = (1:2))
dev.off()


#################################

