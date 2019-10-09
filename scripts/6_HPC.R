#!/usr/bin/Rscript

# language: R 
# version this script was developed in: "R version 3.6.0 (2019-04-26)"
# platform this script was developed in: "x86_64-apple-darwin15.6.0"

#NB: this script doesn't run for: NegativeControl, MockCommunity, Protexin and Coliguard

#This script requires the following packages:
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("plyr", repos = "http://cran.us.r-project.org")
install.packages("seqinr", repos = "http://cran.us.r-project.org")

#upload all libraries
library(seqinr)
library(data.table)
library(plyr)

#create empty df with all possible headers (dates) and interate through each sample's df
#containing depths from clustered weighted bins 
#concatenate all to one df > output as: df_total.csv
pig.id.basedir = "/shared/homes/12705859/out_new"
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
    
  df = data.frame(
    pig = character(),
    bin = character(),
    secondary_cluster = character(),
    X170130.01.bam = character(),
    X170131.01.bam = character(),
    X170131.02.bam = character(),
    X170131.03.bam = character(),
    X170201.01.bam = character(),
    X170201.02.bam = character(),
    X170203.01.bam = character(),
    X170206.01.bam = character(),
    X170206.02.bam = character(),
    X170207.01.bam = character(),
    X170207.02.bam = character(),
    X170208.01.bam = character(),
    X170210.01.bam = character(),
    X170214.01.bam = character(),
    X170214.02.bam = character(),
    X170214.03.bam = character(),
    X170216.01.bam = character(),
    X170216.02.bam = character(),
    X170217.01.bam = character(),
    X170221.01.bam = character(),
    X170224.01.bam = character(),
    X170228.01.bam = character(),
    X170303.01.bam = character(),
    X170306.01.bam = character(),
    X170307.01.bam = character(),
    X170308.01.bam = character(),
    X170309.01.bam = character(),
    X170310.01.bam = character()
  )
  #if statement to exclude NegativeControl, MockCommunity, Protexin and Coliguard
  if (grepl("[[:digit:]]", pig.id) == TRUE) {
    old <- paste(pig.id.dir, "clustered_wa_bins.csv", sep="/")
    df2 <- read.table(file = old, header = TRUE, sep = ",", row.names = NULL)
    
    new_df <- rbind.fill(df,df2)
    
    fwrite(
      x = new_df,
      file = file.path(pig.id.dir, "new_headers_clustered_wa_bins.csv"),
      row.names=FALSE
    )
  }
}


pig.id.basedir = "/shared/homes/12705859/out_new"
df_total = data.frame()
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  #if statement to exclude NegativeControl, MockCommunity, Protexin and Coliguard
  if (grepl("[[:digit:]]", pig.id) == TRUE) {
    old <- paste(pig.id.dir, "new_clustered_wa_bins.csv", sep="/")
    df2 <- read.table(file = old, header = TRUE, sep = ",", row.names = NULL)
    
    df3 <- data.frame(df2)
    df_total <- rbind(df_total,df3)
    
    fwrite(
      x = df_total,
      file = file.path("/shared/homes/12705859/merged_all_clustered_wa_bins.csv"),
      row.names=FALSE
      )
  }
}

