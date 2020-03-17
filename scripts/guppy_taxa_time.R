
library(wordspace)
library(pheatmap)
library(taxize)
library(ggplot2)


# piggies: 


simplified <- complete %>%
  select(sample_type, guppied_date,branch_width,var_explained,component,PC_position,taxa_simple)

# simplified comes from guppy_XML_process.R

# taking all along except the guppy run with all the samples from all time points 
simplified2 <- simplified %>%
  filter(!guppied_date == "all")  
# now the only piggies are from the DF_piggies_time guppy runs
simplified2$sample_type <- gsub("piggies","DF_piggies_time",simplified2$sample_type)

# branch_width * var_explained = importance
simplified2 <- simplified2 %>%
  mutate(importance = as.numeric(branch_width)*var_explained) %>%
  filter(sample_type=="DF_piggies_time"|sample_type=="groupA"|sample_type=="groupB") %>% # eventual selection of only one guppy run 
  select(guppied_date,taxa_simple,importance) 

both <- simplified2 %>% 
  group_by(taxa_simple,guppied_date) %>%
  summarize(Sum_importance = sum(importance))

# addtop n by date
both <- both %>% 
  group_by(guppied_date) %>%
  top_n(20)

unique(both$taxa_simple) 

both_wide <- both %>%
  pivot_wider(names_from = guppied_date, values_from = Sum_importance) %>%
  select(taxa_simple,Ja31,Fe7,Fe14,Fe21,Fe28,Ma3) %>%
  mutate_all(~replace(., is.na(.), 0))


both_wide <- as.data.frame(both_wide)

##############################

# put order


mytaxa <- both_wide[,1]
mytaxa <- gsub("_"," ",mytaxa)

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
  arrange(., X2,X3,X4,X5,X6,X7,X8,X9) %>%
  select(taxa_simple,Ja31,Fe7,Fe14,Fe21,Fe28,Ma3)
head(out5)


##############################
##############################


# to matrix conversion

rownames(out5) <- out5[,1]
out5[,1] <- NULL
#out_m <- as.matrix(out5[rowSums(out5)>5,])


# ##norm by rows 
# out_m <- as.matrix(out5[rowSums(out5)>5,])
# out_m <- normalize.rows(out_m, method = "euclidean", 
#                         tol = 1e-6, inplace = TRUE)
# pdf("guppy_time_heatmap_importance_normbyrow.pdf")
# pheatmap(out_m, display_numbers = T,
#          #main="Mock community contamination",
#          cluster_rows = F, cluster_cols = F, fontsize_number = 8,
#          fontsize_row = 8)
# dev.off()

##############################

##norm by columns 
out_m <- as.matrix(out5[rowSums(out5)>2,])
out_m <- normalize.cols(out_m, method = "euclidean", 
                        tol = 1e-6, inplace = TRUE)
pdf("guppy_time.pdf")
pheatmap(out_m, display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
dev.off()

# from this hetamap, taxa explaining differnece the most are: 
# 

summary

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


# positive controls : 


simplified <- complete %>%
  select(sample_type, guppied_date,branch_width,var_explained,component,PC_position,taxa_simple)

# simplified comes from guppy_XML_process.R

# taking all along except the guppy run with all the samples from all time points 
simplified2 <- simplified %>%
  filter(guppied_date == "controls")  

# branch_width * var_explained = importance
simplified2 <- simplified2 %>%
  mutate(importance = as.numeric(branch_width)*var_explained) %>%
  select(guppied_date,taxa_simple,importance) 

both <- simplified2 %>% 
  group_by(taxa_simple,guppied_date) %>%
  summarize(Sum_importance = sum(importance))

unique(both$taxa_simple) # 61 taxa! 

both_wide <- both %>%
  pivot_wider(names_from = guppied_date, values_from = Sum_importance) %>%
  mutate_all(~replace(., is.na(.), 0))


both_wide <- as.data.frame(both_wide)

##############################

# put order


mytaxa <- both_wide[,1]
mytaxa <- gsub("_"," ",mytaxa)

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
  select(taxa_simple,controls)
head(out5)


##############################
##############################


# to matrix conversion

rownames(out5) <- out5[,1]
out5[,1] <- NULL


##norm by rows 
out_m <- as.matrix(out5)
out_m <- normalize.cols(out_m, method = "euclidean", 
                        tol = 1e-6, inplace = TRUE)
pdf("guppy_time_heatmap_importance_poscontrols.pdf")
pheatmap(out_m, display_numbers = T,
         #main="Mock community contamination",
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
dev.off()

##############################

##norm by columns 
out_m <- as.matrix(out5[rowSums(out5)>5,])
out_m <- normalize.cols(out_m, method = "euclidean", 
                        tol = 1e-6, inplace = TRUE)
pdf("guppy_time_heatmap_importance_normbycol.pdf")
pheatmap(out_m, display_numbers = T,
         #main="Mock community contamination",
         cluster_rows = F, cluster_cols = F, fontsize_number = 8,
         fontsize_row = 8)
dev.off()


##############################



variance <- both %>%
  group_by(taxa_simple) %>%
  summarise(mean=mean(Sum_importance), sd=sd(Sum_importance)) %>%
  arrange(desc(sd))

View(variance)
# take these first n and look with guppy fat if/how they differ by Cohort

variance %>%
  ggplot(.) +
  geom_bar( aes(x=taxa_simple, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=taxa_simple, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  theme(axis.text.x=element_text(angle=90,hjust=1))



