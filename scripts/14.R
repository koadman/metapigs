# 14.R                                        #
# prepares input file to run on fastspar      #
# input: unnormalized data                    #
# 1. subsets to cohorts and dates of interest #
# 2. widest                                   #
# 3. aggregate                                #
# 4. NA to zeros                              #
# 5. rownames and column names formatting     #
# 6. floats to integers                       #
# 7. saved input file for fastspar            #
# ps: if run as is, fine. If it's modified    #
# make sure to run on the command line:       #
# tr '\r' '\n' < old_file > new_file          #
# before feeding it to fastspar               # 

# loads input file: tot_counts_dereplicated.csv from script 7.R (unnormalized data)
tot_counts_dereplicated <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/tot_counts_dereplicated.csv")

library(tidyverse)
# subset
Ctrl_neo_0131_0207 <- tot_counts_dereplicated %>% filter(
  date == "2017-01-31" | date == "2017-02-07", 
  cohort == "Neomycin" | cohort == "Control"
)

# make widest, one row per secondary_custer + bin
library(tidyverse)
Ctrl_neo_0131_0207_widest <- Ctrl_neo_0131_0207 %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")

# average of all bins that belong to the same cluster
bbb <- Ctrl_neo_0131_0207_widest
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("secondary_cluster")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)

#replace missing values with zeros
ccc[is.na(ccc)] <- 0

# save count data and bring secondary_cluster column to rownames
ccc2 <- ccc[,-1]
rownames(ccc2) <- ccc[,1]

#round numbers : floats to int 
ccc3 <- ceiling(ccc2)

# rownames to cols
library(data.table)
ccc4 <- setDT(ccc3, keep.rownames = TRUE)[]

# save as dataframe
ccc5 <- as.data.frame(ccc4)

# sample names (colnames) replaced with sequential numbers except for secondary_cluster column ID --> OTU_id
colnames(ccc5) <- paste0(1:ncol(ccc5))
# nb: header is also case sensitive! must say "OTU_id" !
colnames(ccc5)[1] <- "OTU_id"

# to matrix
ccc5 <- as.matrix(ccc5)

typeof(ccc5)
typeof(original_fake_data)

class(ccc5)
class(original_fake_data)
#different 

# save my df into a csv and test whether when you open it its class is the same as original_fake_data
write.table(x = ccc5, file = "~/Desktop/ccc5.txt", row.names = FALSE, quote = FALSE, sep = '\t')

ccc5_opened <- read_csv("~/Desktop/ccc5.txt")
identical(class(original_fake_data), class(ccc5_opened) )
# if true, the saved file is good to go! 
