
#merge df_total.csv with cohorts.xlsx from metapigs/source_data

library(readxl)

#input files
cohorts <- read_excel("~/metapigs/source_data/cohorts.xlsx")
df_total <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/merged_all_clustered_wa_bins.csv")

colnames(df_total)[colnames(df_total)=="X170130.01.bam"] <- "17-01-30.1"
colnames(df_total)[colnames(df_total)=="X170131.01.bam"] <- "17-01-31.1"
colnames(df_total)[colnames(df_total)=="X170131.02.bam"] <- "17-01-31.2"
colnames(df_total)[colnames(df_total)=="X170131.03.bam"] <- "17-01-31.3"
colnames(df_total)[colnames(df_total)=="X170201.01.bam"] <- "17-02-01.1"
colnames(df_total)[colnames(df_total)=="X170201.02.bam"] <- "17-02-01.2"
colnames(df_total)[colnames(df_total)=="X170203.01.bam"] <- "17-02-03.1"
colnames(df_total)[colnames(df_total)=="X170206.01.bam"] <- "17-02-06.1"
colnames(df_total)[colnames(df_total)=="X170206.02.bam"] <- "17-02-06.2"
colnames(df_total)[colnames(df_total)=="X170207.01.bam"] <- "17-02-07.1"
colnames(df_total)[colnames(df_total)=="X170207.02.bam"] <- "17-02-07.2"
colnames(df_total)[colnames(df_total)=="X170208.01.bam"] <- "17-02-08.1"
colnames(df_total)[colnames(df_total)=="X170210.01.bam"] <- "17-02-10.1"
colnames(df_total)[colnames(df_total)=="X170214.01.bam"] <- "17-02-14.1"
colnames(df_total)[colnames(df_total)=="X170214.02.bam"] <- "17-02-14.2"
colnames(df_total)[colnames(df_total)=="X170214.03.bam"] <- "17-02-14.3"
colnames(df_total)[colnames(df_total)=="X170216.01.bam"] <- "17-02-16.1"
colnames(df_total)[colnames(df_total)=="X170216.02.bam"] <- "17-02-16.2"
colnames(df_total)[colnames(df_total)=="X170217.01.bam"] <- "17-02-17.1"
colnames(df_total)[colnames(df_total)=="X170221.01.bam"] <- "17-02-21.1"
colnames(df_total)[colnames(df_total)=="X170207.01.bam"] <- "17-02-7.1"
colnames(df_total)[colnames(df_total)=="X170224.01.bam"] <- "17-02-24.1"
colnames(df_total)[colnames(df_total)=="X170228.01.bam"] <- "17-02-28.1"
colnames(df_total)[colnames(df_total)=="X170303.01.bam"] <- "17-03-3.1"
colnames(df_total)[colnames(df_total)=="X170306.01.bam"] <- "17-03-6.1"
colnames(df_total)[colnames(df_total)=="X170307.01.bam"] <- "17-03-7.1"
colnames(df_total)[colnames(df_total)=="X170308.01.bam"] <- "17-03-8.1"
colnames(df_total)[colnames(df_total)=="X170309.01.bam"] <- "17-03-9.1"
colnames(df_total)[colnames(df_total)=="X170310.01.bam"] <- "17-03-10.1"

library(dplyr)
library(data.table)

#merge cohorts
df_total_2 <- merge.data.frame(df_total, cohorts, by.x="pig", by.y = "Animal ID")

#move Cohort Name column to first position of dataframe
df_total_3 <- df_total_2 %>% 
  select("Cohort Name", everything())
#rename to eliminate space
colnames(df_total_3)[colnames(df_total_3)=="Cohort Name"] <- "cohort"
#fill out column "secondary_cluster" -empty cells with NA
library(dplyr)
df_total_4 <- mutate_all(df_total_3, funs(na_if(.,"")))

#subset to contain only rows that have a secondary_cluster value (not na)
clustered_only<-subset(df_total_4, (!is.na(df_total_4[,4])))

View(clustered_only)

#subset based on number of occurrences of cluster IDs. 
#which means: if the cluster contains >100 bins (i.e. bins appear in >100 samples) then bring forward with this subset (newdata)
x <- 50 # number of occurrences
y <- split(clustered_only,f=clustered_only$secondary_cluster)
z <- y[unlist(lapply(y,nrow))>x]
newdata <- vector()
for( k in z ) {
  newdata <- rbind(newdata,k)
}

#these are to test if the subsetting to x number of occurrences is working fine
#change the secondary_cluster ID to check which ones are broughtforward in the analysis
#sum(clustered_only$secondary_cluster == "1000_1")
#sum(newdata$secondary_cluster == "1000_1")

#what's the number of uniq clusters we have now? 
#library(data.table)
#DT <- data.table(newdata)

#when we keep only bins that occur 100 times (100 times the same cluster) we end up with 24 different clusters
#DT2 <- DT[, .(number_of_distinct_clusters = length(unique(secondary_cluster))), by = secondary_cluster]
#View(DT2)

#make long (one value per row)
library(reshape)
newdata_long <- melt(newdata, id=c("cohort", "pig", "bin", "secondary_cluster"))

#split column "variable" using dot separator
#install.packages("splitstackshape")
library(splitstackshape)
NL2 <- cSplit(newdata_long, "variable", ".")
# rename new columns
colnames(NL2)[colnames(NL2)=="variable_1"] <- "date"
colnames(NL2)[colnames(NL2)=="variable_2"] <- "replicate"

#reformat to be recognized as dates
#NL2$date <- as.Date(NL2$date, format = "%y-%m-%d")

#keep only rows that have values (exclude NA rows)
all_but_NA <- NL2[complete.cases(NL2), ]
View(all_but_NA)
#return rows that have replicate==1, ==2, ==3, ==4
NLrep1 <- all_but_NA[which(all_but_NA$replicate == 1), ]
NLrep2 <- all_but_NA[which(all_but_NA$replicate == 2), ]
NLrep3 <- all_but_NA[which(all_but_NA$replicate == 3), ]
NLrep4 <- all_but_NA[which(all_but_NA$replicate == 4), ]

#concatenate the above dataframes
total <- rbind(NLrep1, NLrep2, NLrep3, NLrep4)
View(total)
#take mean out of the replicates (works right)
library(dplyr)
total_dereplicated <- total %>% group_by(pig, bin, secondary_cluster, cohort, date) %>% 
  summarize(value.mean = mean(value))
View(total_dereplicated)

#write out to play with it interactively
fwrite(x = total_dereplicated, file = "~/Desktop/bins_clustering_parsing_dataframes/total_dereplicated.csv")

#make sure column date is seen as Date
# (it matters when geom_smooth)
total_dereplicated$date <- as.Date(total_dereplicated$date, format='%y-%m-%d')
class(total_dereplicated$date)


# for uniq secondary_cluster ID, make a single plot: 
# all cohorts in one plot
# works fine, but doesn't look good
library(ggplot2)
# Design a function
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=value.mean))+
    geom_point(aes(colour = factor(cohort)), size = 0.2)+
    #geom_line to connect the dots by pig
    geom_line(aes(group = pig))+
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(total_dereplicated$secondary_cluster), gg_fun, dt = total_dereplicated)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/plots_by_secondary_cluster.pdf")
plot_list
dev.off()

############################################################################################################

# PLOTS: 

#plot multiple plots (1 per cohort) in each page (1 page : 1 cluster) 
# IT WORKS! 
# works as well with addition of smooth
#loess is the stats method chosen if we first choose automatic
# For method = "auto" the smoothing method is chosen based on the size of the largest group (across all panels).
#better to set "loess" than "auto" in case it forces some clusters
#to be done with a different method
library(ggplot2)
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=value.mean, colour = cohort))+
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
plot_list <- lapply(unique(total_dereplicated$secondary_cluster), gg_fun, dt = total_dereplicated)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/plots_by_secondary_cluster.pdf")
plot_list
dev.off()


############################################################################################################

#ggplot playing with just one cluster:
#subset whole dataset to one cluster only
subset <- total_dereplicated[which(total_dereplicated$secondary_cluster == "34_1"),names(total_dereplicated) %in% 
                       c("pig","bin","cohort", "date", "value.mean")]
#subset one cluster one cohort
subset2 <- subset[which(subset$cohort == "Control"),names(subset) %in% 
                               c("pig","bin", "date", "value.mean")]
View(subset2)
#simple and works
ggplot(subset2, aes(date, value.mean)) + 
  geom_point() 

#now do smooth: 
gg2 <- ggplot(subset, aes(date, value.mean)) + 
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
     ~date, value = "value.mean")
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
  #secondary_cluster == "34_1"
)
View(Ctrl_neo_0131_0207)

#make widest, one row per secondary_custer + bin
library(tidyverse)
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value.mean")
View(Ctrl_neo_0131_0207_widest)


#join first two columns "bin + secondary_cluster" to create genome ID
Ctrl_neo_0131_0207_widest$sec_clu_bin <- paste(Ctrl_neo_0131_0207_widest$secondary_cluster, Ctrl_neo_0131_0207_widest$bin, sep="_")
#remove column bin and secondary cluster
Ctrl_neo_0131_0207_widest$bin <- NULL
Ctrl_neo_0131_0207_widest$secondary_cluster <- NULL
#bring the sec_clu_bin column to first position
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207_widest[,c(ncol(Ctrl_neo_0131_0207_widest),1:(ncol(Ctrl_neo_0131_0207_widest)-1))]
View(Ctrl_neo_0131_0207_widest)

#replace NAs with zeros (yet to be determined if it's a good idea in our case)
Ctrl_neo_0131_0207_widest[is.na(Ctrl_neo_0131_0207_widest)] <- 0

fwrite(Ctrl_neo_0131_0207_widest, file = "~/Desktop/bins_clustering_parsing_dataframes/Ctrl_neo_0131_0207_widest.csv")
