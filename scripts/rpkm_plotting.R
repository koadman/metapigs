#to be run before this: 
#./metapigs/scripts/add_extraction_plates_data.py 
#/Users/12705859/metapigs/source_data/pigs_samples_IDs_for_NCBI.xlsx 
#/Users/12705859/metapigs/source_data/DNA_plates.xlsx 
#/Users/12705859/metapigs/source_data/cohorts.xlsx > new_table

#input files for the script below: 
#new_table (created as above)
#rpkm_results (created with rpkm.py)

#rememebr to import the datasets this way: 
#new_table: from text (base)
#rpkm_results: from text (readr) choosing Tab as separator (if imported as (base) the number of entries for some reason will be halved )

setwd("~/Desktop/rpkm_results_dir")

#rememebr to import the datasets this way: 
#new_table: from text (base)
#rpkm_results: from text (readr) choosing Tab as separator (if imported as (base) the number of entries for some reason will be halved )
#or import datasets using code: 
library(readr)
new_table <- read.delim("~/Desktop/rpkm_results_dir/new_table")
View(new_table)
rpkm_results <- read_delim("rpkm_results", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
View(rpkm_results)

#transform to dataframes
df1=data.frame(new_table)
df2=data.frame(rpkm_results)

#give header to rpkm_results. 
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
View(df4)

#first replicate the column (to make sure the changing of the names in that column is working right)
df4$X.collection_date2 = df4$X.collection_date
#then change the (created) column attributes
levels(df4$X.collection_date2) <- c("01/30", "01/31", "02/01", "02/03", "02/06", "02/07", "02/08", "02/10", "02/14", "02/16", "02/17", "02/21", "02/24", "02/28", "03/03", "03/06", "03/07", "03/08", "03/09", "03/10", "08/14", "18/01/24", "NaT")

#subset of useful columns (good for now, for debugging)
df9 <- df4[, c(1, 2, 3, 22, 31, 35)]

#subset to Neomycin only
df10 <- df9[df9$Cohort == "Neomycin", ]

#plotting all the piglets, with a break. loop works but it s not using the right data. 
pig <- levels(df4$isolation_source)
par(mfrow = c(2,1))
cnt <- 1

for (i in df10){
  plot(df10$X.collection_date2, df10$rpkm,
       main = paste(i),
       xlab = "time",
       ylab = "rpkm")
  cnt <- cnt + 1
  if (cnt > 2) {
    break
  }
}


#subset to Neomycin Cohort only
df6 <- subset(df4, Cohort == "Neomycin")
df7 <- subset(df6, gene_hit == "aph(6)-Id_5_18676889_1")

#plotting for all piglets inside subset
par(mfrow = c(3,3))
cnt <- 1

for (i in df7$isolation_source){
  plot(df7$X.collection_date2, df7$rpkm,
       main = paste(i),
       xlab = "time",
       ylab = "rpkm")
  cnt <- cnt + 1
  if (cnt > 7) {
    break
  }
}

plot(df4$X.collection_date2, df4$rpkm)


















  
#plot all genes
plot(df4$Cohort, df4$rpkm, main="aminoglycoside_genes_rpkm", ylab="rpkm")

#to subset it to one specific gene:
df5 <- subset(df4, gene_hit == 'neomycinresistanceprotein"/protein_id="AAA88361.1"', select = c("rpkm", "Cohort", "X.collection_date"))

#plot one gene, setting y (rpkm) limit 
plot(df5$Cohort, df5$rpkm, main="neomycin_res_AAA88361.1_rpkm", ylab="rpkm", ylim=c(0,30))

#to look at which ones of the rpkm in df5 are really high
df5 %>% filter(rpkm > 400)

#order Cohorts on the x axis
df4$Cohort<-factor(df4$Cohort, levels=c("Mothers", "Control", "ColiGuard", "D-scour", "Neomycin", "Neomycin+ColiGuard", "Neomycin+D-scour"))
df5$Cohort<-factor(df5$Cohort, levels=c("Mothers", "Control", "ColiGuard", "D-scour", "Neomycin", "Neomycin+ColiGuard", "Neomycin+D-scour"))

#look at how many hits depending on time 
#filter out the Neomycin cohort into new dataframe
#not working. collection dates are all NA
df7 <- subset(df4, Cohort == 'Neomycin', select = c("rpkm", "Cohort", "X.collection_date"))


#order based on date (otherwise not ordered logically)
df7$X.collection_date<-factor(df7$X.collection_date, levels=c("01/30", "01/31", "02/01", "02/03", "02/06", "02/07", "02/08", "02/10", "02/14", "02/16", "02/17", "02/21", "02/24", "02/28", "03/03", "03/06", "03/07", "03/08", "03/09", "03/10", "08/14", "18/01/24", "NaT" ))

plot(df7$X.collection_date, df7$rpkm, main="aminoglycoside_genes_rpkm_time")



