
# 7.R script                    #
# Further cleaning,             #
# parsing, and                  #
# transform depth into counts   #

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



############################################################################################################

# TEST HYPOTHESES: 

# subset to dates: Jan31 and Fe7
# subset to cohorts: control and neomycin

#subset whole dataset to dates of interest
Ctrl_neo_0131_0207 <- total_dereplicated %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
  #secondary_cluster == "34_1"
)
View(Ctrl_neo_0131_0207)

library(reshape)
#from long to wide (sooooo much better than pivot! this is fast and efficient!) # wasted lot of time to get this working
# checked and it's correct
Ctrl_neo_0131_0207_wide <- cast(data = b, pig + bin + secondary_cluster + cohort
     ~date, value = "value")
View(Ctrl_neo_0131_0207_wide)

#but how many NA values? just to know...
sum(is.na(Ctrl_neo_0131_0207_wide))
#198/916 when subset is Ctrl_neo_0131_0207

#calculate delta
Ctrl_neo_0131_0207_wide$diff <- (Ctrl_neo_0131_0207_wide$`2017-01-31` - Ctrl_neo_0131_0207_wide$`2017-02-07`)
View(Ctrl_neo_0131_0207_wide)

install.packages("dplyr")
library(dplyr)
group_by(Ctrl_neo_0131_0207_wide, cohort) %>%
  summarise(
    count = n(),
    mean = mean(diff, na.rm = TRUE),
    sd = sd(diff, na.rm = TRUE)
  )
#visualize
#install.packages("ggpubr")
library("ggpubr")
ggboxplot(Ctrl_neo_0131_0207_wide, x = "cohort", y = "diff", 
          color = "cohort", palette = c("#00AFBB", "#E7B800"),
          ylab = "diff", xlab = "cohort")

#before running an indipendent t test on the two cohort groups 
#you need to check whether the two groups are normally distributed
#do this with the Shapiro-Wilk test: 
#if p-value > 0.05 then the groups are NOT normally distributed
with(Ctrl_neo_0131_0207_wide, shapiro.test(diff[cohort == "Control"])) ## W = 0.79465, p-value = 0.007999
with(Ctrl_neo_0131_0207_wide, shapiro.test(diff[cohort == "Neomycin"])) #W = 0.64225, p-value = 4.171e-05
# p value < 0.05 means that the two groups are NOT normally distributed

#normally you stop here because the Shapiro-Wilk test result is telling you that the two groups are not normally distributed
#if you ignore this and go ahead anyway: 

# F-test to test for homogeneity in variances
res.ftest <- var.test(diff ~ cohort, data = Ctrl_neo_0131_0207_wide)
res.ftest
# returns p-value = 0.3582. p-value = 0.3582 is greater than the significance level alpha = 0.05. 
#In conclusion, there is no significant difference between the variances of the two sets of data. 
#Therefore, we can use the classic t-test which assume equality of the two variances.

# Compute t-test
res <- t.test(diff ~ cohort, data = Ctrl_neo_0131_0207_wide, var.equal = TRUE)
res
#returns: p-value = 0.7745

# printing the p-value
res$p.value
# printing the mean
res$estimate
# printing the confidence interval
res$conf.int

############################################################################################################

# TEST HYPOTHESES using EdgeR and DESeq. you'll need to make the dataframe the widest
# start from dataframe Ctrl_neo_0131_0207 and make it the widest you can, 
# keeping unique bin + secondary_cluster as rows

#subset whole dataset to dates of interest
Ctrl_neo_0131_0207 <- total_dereplicated %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
)

#make widest, one row per secondary_custer + bin
library(tidyverse)
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")
View(Ctrl_neo_0131_0207_widest)


#join first two columns "bin + secondary_cluster" to create genome ID
Ctrl_neo_0131_0207_widest$sec_clu_bin <- paste(Ctrl_neo_0131_0207_widest$secondary_cluster, Ctrl_neo_0131_0207_widest$bin, sep="_")
#remove column bin and secondary cluster
Ctrl_neo_0131_0207_widest$bin <- NULL
Ctrl_neo_0131_0207_widest$secondary_cluster <- NULL
#bring the sec_clu_bin column to first position
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207_widest[,c(ncol(Ctrl_neo_0131_0207_widest),1:(ncol(Ctrl_neo_0131_0207_widest)-1))]

#replace NAs with zeros (yet to be determined if it's a good idea in our case)
Ctrl_neo_0131_0207_widest[is.na(Ctrl_neo_0131_0207_widest)] <- 0
View(Ctrl_neo_0131_0207_widest)

fwrite(Ctrl_neo_0131_0207_widest, file = "~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_widest_counts.csv")
