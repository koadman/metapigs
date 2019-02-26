#to be run before this: 
#./metapigs/scripts/add_extraction_plates_data.py 
#/Users/12705859/metapigs/source_data/pigs_samples_IDs_for_NCBI.xlsx 
#/Users/12705859/metapigs/source_data/DNA_plates.xlsx 
#/Users/12705859/metapigs/source_data/cohorts.xlsx > new_table

#input files for the script below: 
#new_table (created as above)
#rpkm_results (created with rpkm.py)


setwd("~/")

#read the two input files
new_table <- read.delim("~/Desktop/new_table")
View(new_table)

rpkm_results <- read.table("~/Desktop/rpkm_results", header=FALSE)
View(rpkm_results)

#transform to dataframes
df1=data.frame(new_table)
df2=data.frame(rpkm_results)

#give header to rpkm_results
names(df2)[1] <- "plate_well"
names(df2)[2] <- "gene_hit"
names(df2)[3] <- "rpkm"

#install packages and import 
install.packages("sqldf")
library(sqldf)
install.packages("tidyr")
library(tidyr)

#join columns containing plate and well ID into one, without removing the original columns
df3 <- unite(df1, col = "plate_well", DNA_plate:DNA_well, sep = "_", remove = FALSE)

#merge dataframes new_table(NCBI) and rpkm_results into new dataframe, based on column called "plate_well"
df4 <- merge(df2, df3, by="plate_well")

#plot all genes
plot(df4$Cohort, df4$rpkm, main="aminoglycoside_genes_rpkm", ylab="rpkm")

#plot one gene, setting y (rpkm) limit 
plot(df5$Cohort, df5$rpkm, main="neomycin_res_AAA88361.1_rpkm", ylab="rpkm", ylim=c(0,30))

#to subset it to one specific gene:
df5 <- subset(df4, gene_hit == 'neomycinresistanceprotein"/protein_id="AAA88361.1"', select = c("rpkm", "Cohort", "X.collection_date"))

#to look at which ones of the rpkm in df5 are really high
df5 %>% filter(rpkm > 400)

#order Cohorts on the x axis
df4$Cohort<-factor(df4$Cohort, levels=c("Mothers", "Control", "ColiGuard", "D-scour", "Neomycin", "Neomycin+ColiGuard", "Neomycin+D-scour"))
df5$Cohort<-factor(df5$Cohort, levels=c("Mothers", "Control", "ColiGuard", "D-scour", "Neomycin", "Neomycin+ColiGuard", "Neomycin+D-scour"))

#look at how many hits dpeneding on time 
#filter out the Neomycin cohort into new dataframe
df7 <- subset(df4, Cohort == 'Neomycin', select = c("rpkm", "Cohort", "X.collection_date"))

#order based on date (otherwise not ordered logically)
df7$X.collection_date<-factor(df7$X.collection_date, levels=c("2017-01-30 00:00:00" "2017-01-31 00:00:00" "2017-02-01 00:00:00" "2017-02-03 00:00:00" "2017-02-06 00:00:00"
"2017-02-07 00:00:00" "2017-02-08 00:00:00" "2017-02-10 00:00:00" "2017-02-14 00:00:00" "2017-02-16 00:00:00"
"2017-02-17 00:00:00" "2017-02-21 00:00:00" "2017-02-24 00:00:00" "2017-02-28 00:00:00" "2017-03-03 00:00:00"
"2017-03-06 00:00:00" "2017-03-07 00:00:00" "2017-03-08 00:00:00" "2017-03-09 00:00:00" "2017-03-10 00:00:00"
"2017-08-14 00:00:00" "2018-01-24 00:00:00" "NaT" ))

plot(df7$X.collection_date, df7$rpkm, main="aminoglycoside_genes_rpkm_time")



