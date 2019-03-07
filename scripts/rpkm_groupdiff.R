#!/usr/bin/env Rscript

#
# WARNING: not yet working correctly!
# script to compute group differences in response to neomycin treatment from rpkm values
#

library(tidyr)
library(readr)

new_table <- read.delim("new_table")
rpkm_results <- read_delim("rpkm_results", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

#transform to dataframes
df1=data.frame(new_table)
df2=data.frame(rpkm_results)

#give header to rpkm_results. 
names(df2)[1] <- "plate_well"
names(df2)[2] <- "gene_hit"
names(df2)[3] <- "rpkm"


#join columns containing plate and well ID into one, without removing the original columns
df3 <- unite(df1, col = "plate_well", DNA_plate:DNA_well, sep = "_", remove = FALSE)

#merge dataframes new_table(NCBI) and rpkm_results into new dataframe, based on column called "plate_well"
df4 <- merge(df2, df3, by="plate_well")

#first replicate the column (to make sure the changing of the names in that column is working right)
df4$X.collection_date2 = df4$X.collection_date

df9 <- df4[, c(1, 2, 3, 22, 31, 35)]


ttester <- function(gene) {
  neo <- bobo$delta_rpkm[bobo$Cohort=="Neomycin" & bobo$gene_hit==gene]
  control <- bobo$delta_rpkm[bobo$Cohort=="Control" & bobo$gene_hit==gene]
  ttt <- 1
  if(length(neo)>2&length(control)>2){
    tttest<-t.test(neo, control)
    ttt<-tttest$p.value
  }
  ttt
}


dodo_t0<-df9[df9$X.collection_date2=="2017-01-31 00:00:00",c(2,3,4,5)]
dodo_t1<-df9[df9$X.collection_date2=="2017-02-07 00:00:00",c(2,3,4,5)]
bobo <- merge(dodo_t0,dodo_t1, by=c("isolation_source","gene_hit","Cohort"))
bobo$delta_rpkm<-bobo$rpkm.y-bobo$rpkm.x

ttt<-sapply(unique(bobo$gene_hit),ttester)

summary(ttt)
