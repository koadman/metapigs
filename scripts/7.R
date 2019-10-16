
# 7.R script                    #
# Further cleaning,             #
# parsing, and                  #
# transform depth into counts   #
# playing with plots            #


library(data.table)
merged_all_clustered_wa_bins_with_cohorts <- read_csv("~/Desktop/bins_clustering_parsing_dataframes/merged_all_clustered_wa_bins_with_cohorts_all.csv")

#subset to contain only rows that have a secondary_cluster value (non NA)
aaa<-subset(merged_all_clustered_wa_bins_with_cohorts, (!is.na(merged_all_clustered_wa_bins_with_cohorts[,5])))
#what's the number of uniq clusters we have now? 
length(unique(aaa$secondary_cluster))

#get the column names before entering function
original_colnames <- colnames(aaa)

# transform depth values into counts
# where {counts = depth*binLen/300} (300 is the read pair length)
A <- function(x) x * aaa[,2] / 300
counts <- data.frame(aaa[1:5], apply(aaa[,6:ncol(aaa)],2, A) )

#return the original colnames to new dataframe
colnames(counts) <- original_colnames

#now we don't need binLen anymore we can remove this column
counts <- counts[,-2]

#make long (one value per row)
library(reshape)
counts_long <- melt(counts, id=c("cohort", "pig", "bin", "secondary_cluster"))

#split column "variable" using dot separator
#install.packages("splitstackshape")
library(splitstackshape)
NL2 <- cSplit(counts_long, "variable", ".")
# rename new columns
colnames(NL2)[colnames(NL2)=="variable_1"] <- "date"
colnames(NL2)[colnames(NL2)=="variable_2"] <- "replicate"
#reformat to be recognized as dates
NL2$date <- as.Date(NL2$date, format = "%y-%m-%d")

#keep only rows that have values (exclude NA rows)
all_but_NA <- NL2[complete.cases(NL2), ]

#return rows that have replicate==1, ==2, ==3, ==4
NLrep1 <- all_but_NA[which(all_but_NA$replicate == 1), ]
NLrep2 <- all_but_NA[which(all_but_NA$replicate == 2), ]
NLrep3 <- all_but_NA[which(all_but_NA$replicate == 3), ]
NLrep4 <- all_but_NA[which(all_but_NA$replicate == 4), ]

#concatenate the above dataframes
total <- rbind(NLrep1, NLrep2, NLrep3, NLrep4)

#take mean out of the replicates (works fine, example below)
total_dereplicated <- aggregate(value~cohort+pig+bin+secondary_cluster+date,total,mean)
View(total_dereplicated)

# checking if mean of replicates is working fine
#a <- filter(NLrep1, pig == "29951", bin == "bins.1.fa", date == "2017-01-31")
#a
#b <- filter(NLrep2, pig == "29951", bin == "bins.1.fa", date == "2017-01-31")
#b
#c <- filter(total_dereplicated, pig == "29951", bin == "bins.1.fa", date == "2017-01-31")
#c
# yes, c is the mean of a and b

#write out to play with it 
fwrite(x = total_dereplicated, file = "~/Desktop/bins_clustering_parsing_dataframes/total_dereplicated.csv")

############################################################################################################

# Playing with plots: 
# Let the subsetting and stats begin!

#make sure column date is seen as Date
# (it matters when geom_smooth)
total_dereplicated$date <- as.Date(total_dereplicated$date, format='%y-%m-%d')
class(total_dereplicated$date)

#how many unique clusters is there? so I know how many to expect from the function below? 
length(unique(total_dereplicated$secondary_cluster))

#if secondary_cluster ID occurs more than 100 times, keep rows
library(dplyr)
mostPopular <- total_dereplicated %>%
  group_by(secondary_cluster) %>%
  filter(n() > 500)
View(mostPopular)
nrow(mostPopular)

length(unique(total_dereplicated$secondary_cluster))
# there is 4274 unique secondary clusters in total
length(unique(mostPopular$secondary_cluster))
# there is 52 unique most popular secondary clusters

############################################################################################################

# PLOT 1: 

# for uniq secondary_cluster ID, make a single plot: 
# all cohorts in one plot
# works fine, but doesn't look good
library(ggplot2)
# Design a function
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=value))+
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
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=value, colour = cohort))+
    annotate("rect", 
             xmin = as.Date('2017-02-01'), 
             ymin = -Inf,
             xmax = as.Date('2017-02-05'), 
             ymax = Inf,
             fill = "palegreen") +
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

############################################################################################################

# GGPLOT playing with just one cluster and one cohort:

#subset whole dataset to one cluster only
clu1_1 <- subset(total_dereplicated, secondary_cluster == "1_1", select = c("pig","bin","cohort", "date", "value"))

# test: if same, it's working
#NROW(clu1_1)
#NROW(which(total_dereplicated$secondary_cluster == "1_1"))

#subset further to one cohort
clu1_1_Ctrl <- subset(clu1_1, cohort == "Control", select = c("pig","bin", "date", "value"))

# test: if same, it's working
#NROW(clu1_1_Ctrl)
#NROW(which(clu1_1$cohort == "Control"))

#simple and works
ggplot(clu1_1_Ctrl, aes(date, value)) + 
  geom_point() 

#now do smooth: 
gg2 <- ggplot(clu1_1_Ctrl, aes(date, value)) + 
  annotate("rect", 
           xmin = as.Date('2017-02-01'), 
           ymin = -Inf,
           xmax = as.Date('2017-02-05'), 
           ymax = Inf,
           fill = "palegreen") +
  geom_point(stat='identity', position='identity',size=0.3) 
gg3 <- gg2 + geom_smooth(method='loess', 
              formula = y ~ splines::bs(x, 3))
                        
############################################################################################################

#TO DOs improvements: 
#make rectangles to cover the treatment periods. 
# in plots where one plot : one cohort, you'll
# need to have different rectangles in the same pdf page
#reorder cohorts plots order of appeareance 
#leave mothers out?

#improvement ideas:
#possibly: join column dates together that have +-1 day difference
#interactive plots? 
