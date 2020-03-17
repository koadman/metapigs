library(ggplot2)
library(wordspace)
library(pheatmap)
library(taxize)

setwd("/Users/12705859/Desktop/metapigs_base/phylosift/guppy")
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"

# visualizing guppy fat output 

# load metadata 
mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

# formatting metadata column names 
mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == '*sample_name'] <- 'sample_name'

mdat <- mdat %>%
  dplyr::select(isolation_source,collection_date,Cohort,DNA_plate,DNA_well,PigPen,sample_name)
head(mdat)

# load guppy fat output df

guppyfat_simplified <- read_csv("guppyfat_simplified.df",col_names = TRUE)


# merge guppy fat output df with metadata 

NROW(guppyfat_simplified)
NROW(mdat)

df <- right_join(guppyfat_simplified,mdat)

NROW(df)
head(guppyfat_simplified)

df <- df %>% 
  select(name,branch_length,branch_width,isolation_source,collection_date,Cohort,sample_name)
View(df)
which(is.na(df$name))

##########################################################################################
##########################################################################################

# POSITIVE CONTROLS 

# number of unique taxa identified in whole dataset from (guppy fat) reads placements onto tree
x <- df %>%
  filter(Cohort=="MockCommunity"|Cohort=="PosControl_D-scour"|Cohort=="PosControl_ColiGuard")%>%
  group_by(sample_name) %>% 
  top_n(n = 15, wt = branch_width) %>% # sensible number to not overcrowd heatmap
  select(name,sample_name)

x <- x %>% 
  group_by(sample_name,name) %>% 
  tally() %>%
  pivot_wider(names_from = sample_name, values_from = n) %>%
  mutate_all(~replace(., is.na(.), 0))

x <- as.data.frame(x)

# to matrix conversion
rownames(x) <- x[,1]
x[,1] <- NULL

##norm by cols
x_m <- as.matrix(x)
#x_m <- as.matrix(x_m[rowSums(x_m)>1,])
x_m <- normalize.cols(x_m, method = "euclidean", 
                        tol = 1e-6, inplace = TRUE)
pdf("guppyfat_positive_controls.pdf")
pheatmap(x_m, display_numbers = T,
         #main="Mock community contamination",
         cluster_rows = F, cluster_cols = F, fontsize_number = 3,
         fontsize_row = 3)
dev.off()





##########################################################################################
##########################################################################################

# PIGGIES TIME  -  COHORTS

# number of unique taxa identified in whole dataset from (guppy fat) reads placements onto tree
x <- df %>%
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  top_n(n = 600, wt = branch_width)

# splitting into multiple dataframes (by sample_type name)

multi_x <- split( x , f = x$Cohort )

pdf("guppyfat.pdf")
for (single.x in multi_x) {
  
  which_cohort <- single.x$Cohort
  
  newdf <- single.x %>%
  filter(collection_date=="2017-01-31"|collection_date=="2017-02-07"|collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|collection_date=="2017-02-28"|collection_date=="2017-03-03") %>%
  select(name,collection_date,branch_width)
  
  xxx <- newdf %>% 
    group_by(collection_date,name) %>% 
    summarize(Sum_branch_width = sum(branch_width)) %>%
    pivot_wider(names_from = collection_date, values_from = Sum_branch_width) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  newdf <- as.data.frame(xxx)
  
  # to matrix conversion
  rownames(newdf) <- newdf[,1]
  newdf[,1] <- NULL
  
  ##norm by cols
  x_m <- as.matrix(newdf)
  x_m <- round(normalize.cols(x_m, method = "euclidean", 
                        tol = 1e-6, inplace = TRUE),4)
  
  pheatmap(x_m, display_numbers = T,
           main=as.character(which_cohort)[1],
           cluster_rows = F, cluster_cols = F, fontsize_number = 5,
           fontsize_row = 5)
  
}
dev.off()


##########################################################################################
##########################################################################################

# PIGGIES TIME  - ALL COHORTS

# number of unique taxa identified in whole dataset from (guppy fat) reads placements onto tree
x <- df %>%
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  top_n(n = 600, wt = branch_width)

newdf <- x %>%
  filter(collection_date=="2017-01-31"|collection_date=="2017-02-07"|collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|collection_date=="2017-02-28"|collection_date=="2017-03-03") %>%
  select(name,collection_date,branch_width) %>% 
  group_by(collection_date,name) %>% 
  summarize(Sum_branch_width = sum(branch_width)) %>%
  pivot_wider(names_from = collection_date, values_from = Sum_branch_width) %>%
  mutate_all(~replace(., is.na(.), 0))
  
newdf <- as.data.frame(xxx)
  
# to matrix conversion
rownames(newdf) <- newdf[,1]
newdf[,1] <- NULL
  
##norm by cols
x_m <- as.matrix(newdf)
x_m <- round(normalize.cols(x_m, method = "euclidean", 
                            tol = 1e-6, inplace = TRUE),4)
  

pdf("guppyfat_time.pdf")
pheatmap(x_m, display_numbers = T,
         main=as.character(which_cohort)[1],
         cluster_rows = F, cluster_cols = F, fontsize_number = 5,
         fontsize_row = 5)
dev.off()


##########################################################################################
##########################################################################################

# abundance of probiotic strains in Cohorts 

# subset to cohorts
x <- df %>%
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") 

# reordering
x$Cohort <- factor(x$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

presence <- x[grep("BIFIDOBACTERIACAE", x$name), ]

presence <- x[grep("SALIVARIUS", x$name), ]

presence <- x[grep("LACTOBACILLUS", x$name), ]

unique(x$name)

ggplot(presence, aes(x=collection_date, y=branch_width, group=Cohort, color=Cohort)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "LACTOBACILLUS abundance from branch_width") +
  facet_wrap(~Cohort)


##########################################################################################
##########################################################################################

# PIGGIES TIME  - PHYLA

# number of unique taxa identified in whole dataset from (guppy fat) reads placements onto tree
x <- df %>%
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  filter(collection_date=="2017-01-31"|collection_date=="2017-02-07"|collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|collection_date=="2017-02-28"|collection_date=="2017-03-03") %>%
  select(name,collection_date,branch_width) %>% 
  group_by(collection_date,name) %>% 
  summarize(Sum_branch_width = sum(branch_width)) %>%
  pivot_wider(names_from = collection_date, values_from = Sum_branch_width) %>%
  mutate_all(~replace(., is.na(.), 0))

x <- as.data.frame(x)
head(x)
unique(x$name)
##############################


# put order


x <- x[,1]
x <- gsub("_"," ",x)

uids <- get_uid(mytaxa)

out <- classification(uids)

out2 <- do.call(rbind.data.frame, out)

out2$index <- rownames(out2)

out2 <- cSplit(out2, "index",".")
out2$rank <- NULL
out2$id <- NULL

out3 <- out2 %>% 
  spread(key = index_2, value = name)

out3$index_1 <- NULL
out3 <- as.data.frame(lapply(out3, function(y) gsub(" ", "__", y)))

#remove rows where all NA
out3 <- out3 %>% filter_all(any_vars(!is.na(.)))

ind <- !is.na(out3)
out3$taxa_simple <- tapply(out3[ind], row(out3)[ind], tail, 1)

out3$taxa_simple<-toupper(out3$taxa_simple) 

out4 <- left_join(both_wide,out3,by= "taxa_simple")
NROW(out4)
head(out4)

out5 <- out4 %>%
  arrange(., X2,X3,X4,X5,X6,X7,X8,X9,X10) %>%
  select(taxa_simple,Ja31,Fe7,Fe14,Fe21,Fe28,Ma3)
head(out5)


##############################
##############################


