
# 7.R script                            #
# subsets to cluster containing data    #
# depth to counts                       #
# removes samples replicates            #
# normalizes data                       #
# subsets to x number of clusters       #
# playing with plots                    #

library(data.table)
df <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/merged_all_clustered_wa_bins_with_cohorts.csv", 
                                                      na.strings=c("","NA"),
                                                      check.names = FALSE,
                                                      header = TRUE)
NROW(df)
View(df)
# we remove the secondary_cluster info for the moment, as we need to clean the data
df$secondary_cluster <- NULL

# get the column names before entering function
original_colnames <- colnames(df)

# transform depth values into counts
# where {counts = depth*binLen/300} (300 is the read pair length) #2 is binLen 1:5 is the non-depths data
A <- function(x) x * df[,2] / 300
counts <- data.frame(df[1:4], apply(df[,5:ncol(df)],2, A) )
counts
#return the original colnames to new dataframe
colnames(counts) <- original_colnames
#now we don't need binLen anymore we can remove this column
counts <- counts[,-2]

NROW(counts)

#make long (one value per row)
library(reshape)
counts_long <- melt(counts, id=c("cohort", "pig", "bin"))
NROW(counts_long)

# remove NAs (NAs were automatically created when pivoting large)
counts_long <- na.omit(counts_long)
NROW(counts_long)

# Merge the replicates (some samples are duplicate, same day and pig, so bins are also duplicate)
#split column "variable" using dot separator
#install.packages("splitstackshape")
library(splitstackshape)
NL2 <- cSplit(counts_long, "variable", ".")
NROW(NL2)
# rename new columns
colnames(NL2)[colnames(NL2)=="variable_1"] <- "date"
colnames(NL2)[colnames(NL2)=="variable_2"] <- "replicate"
#reformat to be recognized as dates
NL2$date <- as.Date(NL2$date, format = "%y-%m-%d")

View(NL2)

# checking how many replicates we have
#return rows that have replicate==1, ==2, ==3, ==4
NLrep1 <- NL2[which(NL2$replicate == 1), ]
NLrep2 <- NL2[which(NL2$replicate == 2), ]
NLrep3 <- NL2[which(NL2$replicate == 3), ]
NLrep4 <- NL2[which(NL2$replicate == 4), ]
NROW(NLrep1)
NROW(NLrep2)
NROW(NLrep3)
NROW(NLrep4)
#concatenate the above dataframes
total <- rbind(NLrep1, NLrep2, NLrep3, NLrep4)
NROW(total)

# 1st de-replication (take mean out of replicates): 
NL2$replicate <- NULL
NROW(NL2)
no_reps <- aggregate(NL2$value,by=list(cohort=NL2$cohort,pig=NL2$pig,bin=NL2$bin, date=NL2$date),data=NL2,FUN=mean)
NROW(no_reps)
#replace column name x with "value"
names(no_reps)[names(no_reps) == 'x'] <- 'value'


library(dplyr)
# checking if mean of replicates is working fine
a <- filter(NLrep1, pig == "29951", bin == "bins.1.fa", date == "2017-01-31")
a
b <- filter(NLrep2, pig == "29951", bin == "bins.1.fa", date == "2017-01-31")
b
c <- filter(no_reps, pig == "29951", bin == "bins.1.fa", date == "2017-01-31")
c
# yes, c is the mean of a and b

########################################################

# Important change 20191108: because some piglets could not be sampled on given sampling date
# (or low quality sample was taken, so sample was taken again the next day),
# they were re-sampled the morning after (before treatment was administered)
# for this reason we see the following:

# tM <- "2017-01-30"
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-01-30", 
  replacement = "tM", 
  fixed = TRUE)
# t0 <- "2017-01-31" "2017-02-01" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-01-31", 
  replacement = "t0", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-01", 
  replacement = "t0", 
  fixed = TRUE)

# t1 <- "2017-02-03" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-03", 
  replacement = "t1", 
  fixed = TRUE)

# t2 <- "2017-02-06" "2017-02-07" "2017-02-08"
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-06", 
  replacement = "t2", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-07", 
  replacement = "t2", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-08", 
  replacement = "t2", 
  fixed = TRUE)

# t3 <- "2017-02-10" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-10", 
  replacement = "t3", 
  fixed = TRUE)

# t4 <- "2017-02-14"
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-14", 
  replacement = "t4", 
  fixed = TRUE)

# t4 <- "2017-02-16" "2017-02-17" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-16", 
  replacement = "t4", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-17", 
  replacement = "t4", 
  fixed = TRUE)

# t6 <- "2017-02-21" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-21", 
  replacement = "t6", 
  fixed = TRUE)

# t7 <- "2017-02-24" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-24", 
  replacement = "t7", 
  fixed = TRUE)

# t8 <- "2017-02-28" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-02-28", 
  replacement = "t8", 
  fixed = TRUE)

# t9 <- "2017-03-03" 
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-03-03", 
  replacement = "t9", 
  fixed = TRUE)

# t10 <- "2017-03-06" "2017-03-07" "2017-03-08" "2017-03-09" "2017-03-10"
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-03-06", 
  replacement = "t10", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-03-07", 
  replacement = "t10", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-03-08", 
  replacement = "t10", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-03-09", 
  replacement = "t10", 
  fixed = TRUE)
no_reps[4] <- lapply(
  no_reps[4], 
  gsub, 
  pattern = "2017-03-10", 
  replacement = "t10", 
  fixed = TRUE)

View(no_reps)


# need to substitute cohort names that contain symbols: 
no_reps[1] <- lapply(
  no_reps[1], 
  gsub, 
  pattern = "Neomycin+D-scour", 
  replacement = "NeoD", 
  fixed = TRUE)
no_reps[1] <- lapply(
  no_reps[1], 
  gsub, 
  pattern = "Neomycin+ColiGuard", 
  replacement = "NeoC", 
  fixed = TRUE)
no_reps[1] <- lapply(
  no_reps[1], 
  gsub, 
  pattern = "D-scour", 
  replacement = "Dscour", 
  fixed = TRUE)



# 2nd de-replication: now there should be again replicates because we joined some dates. 
# (take mean out of replicates)
NROW(NL2)
no_reps$sampleID <- paste0(no_reps$date,"_",no_reps$cohort,"_",no_reps$pig,"_",no_reps$bin)
NROW(no_reps) == length(unique(no_reps$sampleID))
NROW(no_reps) # 362863
length(unique(no_reps$sampleID)) # 354050
# here above we can see the amount of replicates (same pig, same date)
no_reps2 <- aggregate(no_reps$value,by=list(cohort=no_reps$cohort,
                                            pig=no_reps$pig,
                                            bin=no_reps$bin,
                                            date=no_reps$date),
                      data=no_reps,FUN=mean)

#replace column name x with "value"
names(no_reps2)[names(no_reps2) == 'x'] <- 'value'

# sanity check: 
# number of unique pigs
length(unique(paste0(no_reps2$pig)))
# 168
# number of bins in total across all pigs
length(unique(paste0(no_reps2$pig,no_reps2$bin)))
# 50787
# number of unique pig,date,bin
length(unique(paste0(no_reps2$pig,no_reps2$date, no_reps2$bin)))
# 354050


# final number of samples per date and cohort : 
az <- no_reps2 %>%
  select(date, pig, cohort)
bz <- unique(az)
cz2 <- data.frame(table(bz$cohort, bz$date))
# subset to non zero containing groups 
dz2 <- cz2 %>% filter(Freq != 0 )
View(dz2)

# Write out to play with it (this is with secondary_clusters assigned only)
fwrite(x = tot_counts_dereplicated, file = "~/Desktop/bins_clustering_parsing_dataframes/tot_counts_dereplicated.csv")
View(tot_counts_dereplicated)

# Write out to play with it (this is everything)
fwrite(x = no_reps2, file = "~/Desktop/bins_clustering_parsing_dataframes/no_reps_all.csv")


########################################################


# normalization (value is divided by the sum of the values per sample)
library(dplyr)
normalized <- tot_counts_dereplicated2 %>%
  group_by(pig, date) %>%
  mutate(norm_counts = value/sum(value))

#remove value column
normalized <- normalized[,-6]


# write out the normalized counts
fwrite(x = normalized, file = "~/Desktop/bins_clustering_parsing_dataframes/normalized.csv")

############################################################################################################

# Playing with plots: 
# Let the subsetting and stats begin!

#make sure column date is seen as Date
# (it matters when geom_smooth)
#normalized$date <- as.Date(normalized$date, format='%y-%m-%d')
#class(normalized$date)

#how many unique clusters is there? so I know how many to expect from the function below? 
length(unique(normalized$secondary_cluster))

#if secondary_cluster ID occurs more than 100 times, keep rows
library(dplyr)
mostPopular <- normalized %>%
  group_by(secondary_cluster) %>%
  filter(n() > 400)
View(mostPopular)
nrow(mostPopular)
#106899 rows left when subsetting (from 162440)

length(unique(normalized$secondary_cluster))
# there is 4465 unique secondary clusters in total
length(unique(mostPopular$secondary_cluster))
# there is 338 unique most popular secondary clusters

############################################################################################################

# PLOT 1: 

# for uniq secondary_cluster ID, make a single plot: 
# all cohorts in one plot
# works fine, but doesn't look good
library(ggplot2)
# Design a function
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=norm_counts))+
    geom_point(aes(colour = factor(cohort)), size = 0.2)+
    #geom_line to connect the dots by pig
    geom_line(aes(group = pig))+
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(mostPopular$secondary_cluster), gg_fun, dt = mostPopular)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/plots_by_secondary_cluster_mostPop.pdf")
plot_list
dev.off()

############################################################################################################

# PLOT 2: 

#plot multiple plots (1 per cohort) in each page (1 page : 1 cluster) 
# IT WORKS! 
# works as well with addition of smooth
#loess is the stats method chosen if we first choose automatic
# For method = "auto" the smoothing method is chosen based on the size of the largest group (across all panels).
#better to set "loess" than "auto" in case it forces some clusters
#to be done with a different method
library(ggplot2)
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=norm_counts, colour = cohort))+
    geom_point(stat='identity', position='identity', aes(colour=cohort),size=0.3) + 
    geom_smooth(method='loess', 
                formula = y ~ splines::bs(x, 3)) +
    facet_wrap(~ cohort, ncol = 2, scales = "free_x") +
    guides(colour = "none") +
    #geom_line to connect the dots by pig
    #geom_line(aes(group = pig))
    theme(text = element_text(size=5), axis.text.x = element_text(angle = 45, hjust=1)) +
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(mostPopular$secondary_cluster), gg_fun, dt = mostPopular)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/plots_by_secondary_cluster_small_mostPop.pdf")
plot_list
dev.off()


mostPopular$secondary_cluster

aaa <- subset(mostPopular, mostPopular$secondary_cluster=="957_1")
bbb <- subset(aaa, aaa$cohort %in% c("Control","Neomycin"))

ggplot(bbb, aes(x=date, y=norm_counts, colour = cohort))+
  geom_point(stat='identity', position='identity', aes(colour=cohort),size=0.3) + 
  geom_smooth(method='loess', 
              formula = y ~ splines::bs(x, 3)) +
  facet_wrap(~ cohort, ncol = 2, scales = "free_x") +
  guides(colour = "none") +
  #geom_line to connect the dots by pig
  #geom_line(aes(group = pig))
  theme(text = element_text(size=5), axis.text.x = element_text(angle = 45, hjust=1))

############################################################################################################

# GGPLOT playing with just one cluster and one cohort:

#subset whole dataset to one cluster only
clu1_1 <- subset(mostPopular, secondary_cluster == "1_1", select = c("pig","bin","cohort", "date", "norm_counts"))

#subset further to one cohort
clu1_1_Ctrl <- subset(clu1_1, cohort == "Control", select = c("pig","bin", "date", "norm_counts"))

#simple and works
ggplot(clu1_1_Ctrl, aes(date, norm_counts)) + 
  geom_point() 

#now do smooth: 
gg2 <- ggplot(clu1_1_Ctrl, aes(date, norm_counts)) + 
  geom_point(stat='identity', position='identity',size=0.3) 
gg3 <- gg2 + geom_smooth(method='loess', 
              formula = y ~ splines::bs(x, 3))
gg3                        

############################################################################################################