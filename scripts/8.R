#!/usr/bin/env Rscript

# 8.R script                                              #
# subsets to cohort and dates of interest,                #
# make widest                                             #
# average of bins that belong to the same cluster         #
# t-tests                                                 #


############################################################################################################

# this assumes the script will be run from the top-level directory of the git repo
data_prefix <- "results/"


library(readr)
library(data.table)
library(dplyr)

# load input file (normalized dataset)
normalized <- read.csv(paste(data_prefix,"normalized.csv",sep=""))
View(normalized)


# EdgeR from new manual:

# subset to dates: Jan31 and Fe7
# subset to cohorts: control and neomycin
Ctrl_neo_0131_0207 <- normalized %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
)

#make widest, one row per secondary_custer + bin
library(tidyverse)
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "norm_counts")

# average of all bins that belong to the same cluster (works!)
bbb <- Ctrl_neo_0131_0207_widest
View(bbb)
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("secondary_cluster")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)

#replace missing values with zeros
ccc[is.na(ccc)] <- 0

View(ccc)

fwrite(ccc, file = paste(data_prefix,"Ctrl_neo_0131_0207_widest.csv",sep=""))


# TEST HYPOTHESES: 

# subset to dates: Jan31 and Fe7
# subset to cohorts: control 
Ctrl_0131_0207 <- normalized %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Control"
)
View(Ctrl_0131_0207)

# subset to cohorts: neomycin
Neo_0131_0207 <- normalized %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin"
)
View(Neo_0131_0207)

# Compare a cohort against itself: 
# can't do the pairedData thing because even though the 
# samples are the same when comparing one cohort against
# itself (t=0 vs t=1) the number of rows (datapoints) 
# is different because of the nature of the method 
# a bin goes down and therefore missing 

# maybe one way to approach this is by making the df wider
# (column 1 is t=0; column 2 is t=1) hereby forcing it 
# to have the same lengths., then you can pivot long 
# again and run the paired t-test
library(reshape)

Ctrl_0131_0207_wide <- cast(data = Ctrl_0131_0207, pig + bin + secondary_cluster
                                ~date, value = "norm_counts")
NROW(Ctrl_0131_0207_wide)
sum(is.na(Ctrl_0131_0207_wide))
View(Ctrl_0131_0207_wide)

Neo_0131_0207_wide <- cast(data = Neo_0131_0207, pig + bin + secondary_cluster
                            ~date, value = "norm_counts")
NROW(Neo_0131_0207_wide)
sum(is.na(Neo_0131_0207_wide))
View(Neo_0131_0207_wide)

#make long again
library(reshape)
Ctrl_0131_0207_long <- melt(Ctrl_0131_0207_wide, id=c("cohort", "pig", "bin", "secondary_cluster"))
Neo_0131_0207_long <- melt(Neo_0131_0207_wide, id=c("cohort", "pig", "bin", "secondary_cluster"))
NROW(Ctrl_0131_0207_long)
NROW(Neo_0131_0207_long)
View(Ctrl_0131_0207_long)
View(Neo_0131_0207_long)

#Nas to zeros? 
Ctrl_0131_0207_long[is.na(Ctrl_0131_0207_long)] <- 0
Neo_0131_0207_long[is.na(Neo_0131_0207_long)] <- 0

# Subset counts data before treatment
Ctrl_0131 <- Ctrl_0131_0207_long %>% filter(
  date == "2017-01-31"
)
NROW(Ctrl_0131)
# Subset counts data after treatment
Ctrl_0207 <- Ctrl_0131_0207_long %>% filter(
  date == "2017-02-07"
)
NROW(Ctrl_0207)

# Plot paired data
library(PairedData)
pd <- paired(Ctrl_0131, Ctrl_0207)
plot(pd, type = "profile") + theme_bw()

# Compute t-test
res <- t.test(Ctrl_0131, Ctrl_0207, paired = TRUE)
res
# p-value = 0.003545, significant difference between before and after 

# Subset counts data before treatment
Neo_0131 <- Neo_0131_0207_long %>% filter(
  date == "2017-01-31"
)
NROW(Neo_0131)
# Subset counts data after treatment
Neo_0207 <- Neo_0131_0207_long %>% filter(
  date == "2017-02-07"
)
NROW(Neo_0207)


# Plot paired data
library(PairedData)
dp <- paired(Neo_0131, Neo_0207)
plot(dp, type = "profile") + theme_bw()

# Compute t-test
res <- t.test(Neo_0131, Neo_0207, paired = TRUE)
res

View(Ctrl_0131)
View(Ctrl_0207)
View(Neo_0131)
View(Neo_0207)



# Plot counts by date and color by date: Control: 
library(ggpubr)
ggboxplot(Ctrl_0131_0207_long, x = "date", y = "value", 
          color = "date", palette = c("#00AFBB", "#E7B800"),
          order = c("2017-01-31", "2017-02-07"),
          ylab = "counts", xlab = "date")
# Plot counts by date and color by date: Neomycin: 
library(ggpubr)
ggboxplot(Neo_0131_0207_long, x = "date", y = "value", 
          color = "date", palette = c("#00AFBB", "#E7B800"),
          order = c("2017-01-31", "2017-02-07"),
          ylab = "counts", xlab = "date")

# Compute summary statistics by date:
# Control: 
library(dplyr)
group_by(Ctrl_0131_0207_long, date) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  )

# Compute summary statistics by date:
# Neomycin: 
library(dplyr)
group_by(Neo_0131_0207_long, date) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  )

# check if the two groups follow a normal distribution
# compute the difference
d1 <- with(Ctrl_0131_0207_long, 
          value[date == "2017-01-31"] - value[date == "2017-02-07"])
d2 <- with(Neo_0131_0207_long, 
          value[date == "2017-01-31"] - value[date == "2017-02-07"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d1)
shapiro.test(d2)
# both are significant. it means the two groups are not normally distributed



################################################


# Let's play with the diffs: 

# NA to zeros 
Ctrl_0131_0207_wide[is.na(Ctrl_0131_0207_wide)] <- 0
View(Ctrl_0131_0207_wide)
Neo_0131_0207_wide[is.na(Neo_0131_0207_wide)] <- 0
View(Neo_0131_0207_wide)

# let's calculate the diff within groups and compare between (ctrl vs neo) groups
#calculate delta
Ctrl_0131_0207_wide$diff <- (Ctrl_0131_0207_wide$`2017-01-31` - Ctrl_0131_0207_wide$`2017-02-07`)
View(Ctrl_0131_0207_wide)
Neo_0131_0207_wide$diff <- (Neo_0131_0207_wide$`2017-01-31` - Neo_0131_0207_wide$`2017-02-07`)
View(Neo_0131_0207_wide)

# rbind the two dfs
Ctrl_neo_0131_0207_wide <- rbind(Ctrl_0131_0207_wide, Neo_0131_0207_wide)
View(Ctrl_neo_0131_0207_wide)

#visualize diff
#install.packages("ggpubr")
library("ggpubr")
ggboxplot(Ctrl_neo_0131_0207_wide, x = "cohort", y = "diff", 
          color = "cohort", palette = c("#00AFBB", "#E7B800"),
          ylab = "diff", xlab = "cohort")

#the stuff above (diffs) is incorrect
#we shouldn not even compare diffs between the groups
#we should instead run an indip t test?



################################################



# you can run an indip t test between two cohorts at the same time point
# to see whether there's a significant amount of bins (clusters)
# differentitally expressed (present) at the same timepoint

# 1. subset to date t=0 including ctrl and neo
# run indip t test

# 2. subset to date t=1 including ctrl and neo
# run indip t test


# 1. 
Ctrl_neo_0131 <- normalized %>% filter(
  date == "2017-01-31", 
  cohort == "Control" | cohort == "Neomycin"
)
View(Ctrl_neo_0131)

# 2. 
Ctrl_neo_0207 <- normalized %>% filter(
  date == "2017-02-07", 
  cohort == "Control" | cohort == "Neomycin"
)
View(Ctrl_neo_0207)

# 1. 
library("ggpubr")
ggboxplot(Ctrl_neo_0131, x = "cohort", y = "norm_counts", 
          color = "cohort", palette = c("#00AFBB", "#E7B800"),
          ylab = "counts", xlab = "cohorts at 0131")

# 2. 
library("ggpubr")
ggboxplot(Ctrl_neo_0207, x = "cohort", y = "norm_counts", 
          color = "cohort", palette = c("#00AFBB", "#E7B800"),
          ylab = "counts", xlab = "cohorts at 0207")

#before running an indipendent t test on the two cohort groups 
#you need to check whether the two groups are normally distributed
#do this with the Shapiro-Wilk test: 
with(Ctrl_neo_0131, shapiro.test(norm_counts[cohort == "Control"])) # W = 0.34547, p-value < 2.2e-16
with(Ctrl_neo_0131, shapiro.test(norm_counts[cohort == "Neomycin"])) # W = 0.32672, p-value < 2.2e-16
with(Ctrl_neo_0207, shapiro.test(norm_counts[cohort == "Control"])) # W = 0.42284, p-value < 2.2e-16
with(Ctrl_neo_0207, shapiro.test(norm_counts[cohort == "Neomycin"])) # W = 0.35545, p-value < 2.2e-16
#if p-value is significant, then the groups are NOT normally distributed
#normally you stop here because the Shapiro-Wilk test result is telling you that the two groups are not normally distributed
#if you ignore this and go ahead anyway: 
# F-test to test for homogeneity in variances
res.ftest <- var.test(value ~ cohort, data = Ctrl_neo_0131)
res.ftest
res.ftest <- var.test(value ~ cohort, data = Ctrl_neo_0207)
res.ftest
# returns p-value < 2.2e-16, which is lower than the significance level alpha = 0.05. 
#In conclusion, there is a significant difference between the variances of the two sets of data. 
#Therefore, we shouldn't use the classic t-test which assume equality of the two variances.
# Compute t-test anyway: 
res0131 <- wilcox.test(value ~ cohort, data = Ctrl_neo_0131,
                   exact = FALSE)
res0131
# W = 6703900, p-value = 1.115e-05 -> counts at 0131 between cohorts are sign different
res0207 <- wilcox.test(value ~ cohort, data = Ctrl_neo_0207,
                       exact = FALSE)
res0207
# W = 5010113, p-value = 0.06387 -> not sign different counts at 0207 between cohorts



################################################



# I now understand the diff between dependent (paired) and indip t-test and where to apply them in this case, 
# but what I need is to run SEPARATE t-tests per each secondary_cluster
# compute the diff
# and test if the difference is significantly different between the two cohorts

# subset to dates: Jan31 and Fe7
# subset to cohorts: control and neo
Ctrl_neo_0131_0207 <- normalized %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Control" | cohort == "Neomycin"
)
View(Ctrl_neo_0131_0207)


# 1. wide
# 2. na to zeros 
# 3. diff 
# 4. for each sec_clu run t test on diffs

library(reshape)

Ctrl_neo_0131_0207_wide <- cast(data = Ctrl_neo_0131_0207, pig + bin + secondary_cluster + cohort
                                ~date, value = "norm_counts")
NROW(Ctrl_neo_0131_0207_wide)
sum(is.na(Ctrl_neo_0131_0207_wide))

# 2. NAs to zeros 
Ctrl_neo_0131_0207_wide[is.na(Ctrl_neo_0131_0207_wide)] <- 0
View(Ctrl_neo_0131_0207_wide)

# 3. diff 
Ctrl_neo_0131_0207_wide$diff <- (Ctrl_neo_0131_0207_wide$`2017-01-31` - Ctrl_neo_0131_0207_wide$`2017-02-07`)
View(Ctrl_neo_0131_0207_wide)

# 4. for each sec_clu run t-test on diffs (by cohort)
lapply(split(Ctrl_neo_0131_0207_wide, Ctrl_neo_0131_0207_wide$secondary_cluster),
       function(d) t.test(Ctrl_neo_0131_0207_wide$diff ~ Ctrl_neo_0131_0207_wide$cohort, data = d))



for (i in unique(Ctrl_neo_0131_0207_wide$secondary_cluster)) {
  print(t.test(Ctrl_neo_0131_0207_wide$diff ~ Ctrl_neo_0131_0207_wide$cohort, group_by=i))
}



# maybe another way would be to make mutliple subsets for each secondary_cluster, then run the t-test for each

# for each unique secondary cluster make a subset dataframe
# for each dataframe run a t-test

datalist = list()
uniq <- unique(unlist(Ctrl_neo_0131_0207_wide$secondary_cluster))
for (i in 1:length(uniq)){
  data_1 <- subset(Ctrl_neo_0131_0207_wide, secondary_cluster == uniq[i], select = c("pig","bin","cohort","diff"))
  datalist[[i]] <- data_1
}

big_data = do.call(rbind, datalist)
View(big_data)




datalist = list()
uniq <- unique(unlist(Ctrl_neo_0131_0207_wide$secondary_cluster))
for (i in 1:length(uniq)){
  data_1 <- subset(Ctrl_neo_0131_0207_wide, secondary_cluster == uniq[i], select = c("pig","bin","cohort","diff","secondary_cluster"))
  dat <- t.test(data_1$diff ~ data_1$cohort)
  dat$i <- i
  datalist[[i]] <- dat
}
big_data = do.call(rbind, datalist)
View(big_data)

#improve the way it stores the t-test outputs

sec_clu_1_1 <- subset(Ctrl_neo_0131_0207_wide, secondary_cluster == "1_1", select = c("pig","bin","cohort","diff","secondary_cluster"))
t.test(sec_clu_1_1$diff ~ sec_clu_1_1$cohort)
