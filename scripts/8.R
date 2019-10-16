
# 8.R script                                              #
# subsetting to cohort and dates of interest,             #
# make wide,                                              #
# t-test                                                  #


############################################################################################################

# load input file (whole dataset)
total_dereplicated <- read_csv("~/Desktop/bins_clustering_parsing_dataframes/total_dereplicated.csv")
View(total_dereplicated)


# TEST HYPOTHESES: 

# subset to dates: Jan31 and Fe7
# subset to cohorts: control and neomycin
Ctrl_neo_0131_0207 <- total_dereplicated %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
)

View(Ctrl_neo_0131_0207)


library(reshape)
#from long to wide (sooooo much better than pivot! this is fast and efficient!) 
# wasted lot of time to get this working
# checked and it's correct
Ctrl_neo_0131_0207_wide <- cast(data = Ctrl_neo_0131_0207, pig + bin + secondary_cluster + cohort
                                ~date, value = "value")

View(Ctrl_neo_0131_0207_wide)

#but how many NA values? just to know...
sum(is.na(Ctrl_neo_0131_0207_wide))
#1094/7562 when subset is Ctrl_neo_0131_0207_wide

#calculate delta
Ctrl_neo_0131_0207_wide$diff <- (Ctrl_neo_0131_0207_wide$`2017-01-31` - Ctrl_neo_0131_0207_wide$`2017-02-07`)
View(Ctrl_neo_0131_0207_wide)

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

# returns p-value < 2.2e-16, which is lower than the significance level alpha = 0.05. 
#In conclusion, there is a significant difference between the variances of the two sets of data. 
#Therefore, we can't use the classic t-test which assume equality of the two variances.

# Compute t-test
res <- t.test(diff ~ cohort, data = Ctrl_neo_0131_0207_wide, var.equal = TRUE)
res
#returns: p-value = 0.04396

# printing the p-value
res$p.value
# printing the mean
res$estimate
# printing the confidence interval
res$conf.int





