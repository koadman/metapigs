

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


#subset based on number of occurrences of cluster IDs. 
#which means: if the cluster contains >100 bins (i.e. bins appear in >100 samples) then bring forward with this subset (newdata)
x <- 100 # number of occurrences
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
install.packages("splitstackshape")
library(splitstackshape)
NL2 <- cSplit(newdata_long, "variable", ".")
# rename new columns
colnames(NL2)[colnames(NL2)=="variable_1"] <- "date"
colnames(NL2)[colnames(NL2)=="variable_2"] <- "replicate"

#reformat to be recognized as dates
NL2$date <- as.Date(NL2$date, format = "%y-%m-%d")
NL2$date <- format(NL2$date, "%y-%m-%d")

#keep only rows that have values (exclude NA rows)
all_but_NA <- NL2[complete.cases(NL2), ]

#return rows that have replicate==1, ==2, ==3, ==4
NLrep1 <- all_but_NA[which(all_but_NA$replicate == 1), ]
NLrep2 <- all_but_NA[which(all_but_NA$replicate == 2), ]
NLrep3 <- all_but_NA[which(all_but_NA$replicate == 3), ]
NLrep4 <- all_but_NA[which(all_but_NA$replicate == 4), ]

#concatenate the above dataframes
total <- rbind(NLrep1, NLrep2, NLrep3, NLrep4)

#take mean out of the replicates (works right)
library(dplyr)
total_dereplicated <- total %>% group_by(pig, bin, secondary_cluster, cohort, date) %>% 
  summarize(value.mean = mean(value))
View(total_dereplicated)

#write out to play with it interactively
fwrite(x = total_dereplicated, file = "~/Desktop/bins_clustering_parsing_dataframes/total_dereplicated.csv")




#let's try to plot one of them: 
#for uniq secondary_cluster ID, make a plot: 
#group_by: cohort
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



#plot multiple plots (1 per cohort) in each page (1 page : 1 cluster) IT WORKS! 
library(ggplot2)
gg_fun <- function(parameter, dt){
  p <- ggplot(dt[dt$secondary_cluster == parameter, ], aes(x=date, y=value.mean, colour = factor(cohort)))+
    geom_point() + facet_wrap(~ cohort, ncol = 2, scales = "free") +
    guides(colour = "none") +
    #geom_line to connect the dots by pig
    geom_line(aes(group = pig)) + 
    theme(text = element_text(size=5), axis.text.x = element_text(angle = 45, hjust=1)) +
    ggtitle(parameter)
  return(p)
}
plot_list <- lapply(unique(total_dereplicated$secondary_cluster), gg_fun, dt = total_dereplicated)
#save all plots into one pdf
pdf("~/Desktop/bins_clustering_parsing_dataframes/plots_by_secondary_cluster.pdf")
plot_list
dev.off()




#TO DOs improvements: 
# y axis limit make it depend on secondary_cluster max value per secondary_cluster
#reorder cohorts plots appeareance 
#leave mothers out










#improvement: 
#possibly: join column dates together that have +-1 day difference
#interactive plots? 
