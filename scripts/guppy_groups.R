
######################################################################################################


# 0 # set working directory & load libs

setwd("/Users/12705859/Desktop/metapigs_base/phylosift/guppy")

library(vcd)
library(summarytools)
library(readr)
library(splitstackshape)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)

# 1 # prepares metadata files to feed guppy_epca.sh:
# one group = one time point : 1 breed : max 2 days diff in bdays

# load metadata 
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"
mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

# formatting metadata column names 
mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == '*sample_name'] <- 'sample_name'


# load breed and bday data 
details <- read_excel(paste0(basedir, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'isolation_source'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'stig'
colnames(details)[colnames(details) == '...8'] <- 'breed'
details$isolation_source <- gsub("G","",details$isolation_source)
details$isolation_source <- gsub("T","",details$isolation_source)

details <- details %>%
  dplyr::select(isolation_source,BIRTH_DAY,breed)


###########################################################################################
###########################################################################################


# POSITIVE CONTROLS 

# filter to keep only pos controls
mdat_sel <- mdat %>% 
  filter(Cohort=="MockCommunity"|Cohort=="PosControl_D-scour"|Cohort=="PosControl_ColiGuard") %>% 
  dplyr::select(isolation_source,DNA_plate,DNA_well)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

mdat_sel$isolation_source
# we have 9 Mock , 8 Protexin, 8 ColiGuard = 25 in total 

mdat_sel <- as.character(mdat_sel$ID)

writeLines(unlist(mdat_sel), "guppy_input/pos_controls_sel.txt", sep = " ")
# contains 25 sample IDs in total 


###########################################################################################
###########################################################################################


# ALL PIGGIES 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

df1 <- setDT(mdat_sel)[, .(Freq = .N), by = .(Cohort)]
df1[order(df1$Cohort)]
# we have 142 ctrl, 133 d-scour, 131 coliguard, 150 neo, 139 neo+d-scour, 130 neo+coliguard

mdat_sel <- as.character(mdat_sel$ID)

NROW(mdat_sel)
# contains 825 sample IDs in total 

writeLines(unlist(mdat_sel), "guppy_input/piggies_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


# ALL PIGGIES ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_Ja31 <- mdat_sel %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

mdat_Fe7 <- mdat_sel %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

mdat_Fe14 <- mdat_sel %>% 
  filter(collection_date == "2017-02-14")

mdat_Fe21 <- mdat_sel %>% 
  filter(collection_date == "2017-02-21")

mdat_Fe28 <- mdat_sel %>% 
  filter(collection_date == "2017-02-28")

mdat_Ma3 <- mdat_sel %>% 
  filter(collection_date == "2017-03-03")

a <- as.character(mdat_Ja31$ID)
b <- as.character(mdat_Fe7$ID)
c <- as.character(mdat_Fe14$ID)
d <- as.character(mdat_Fe21$ID)
e <- as.character(mdat_Fe28$ID)
f <- as.character(mdat_Ma3$ID)

writeLines(unlist(a), "guppy_input/piggies_Ja31_sel.txt", sep = " ")
writeLines(unlist(b), "guppy_input/piggies_Fe7_sel.txt", sep = " ")
writeLines(unlist(c), "guppy_input/piggies_Fe14_sel.txt", sep = " ")
writeLines(unlist(d), "guppy_input/piggies_Fe21_sel.txt", sep = " ")
writeLines(unlist(e), "guppy_input/piggies_Fe28_sel.txt", sep = " ")
writeLines(unlist(f), "guppy_input/piggies_Ma3_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


# SUBSET PIGGIES (two breeds, 3 days bday diff within each breed); groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# merge breed and bday details to metadata
mdat_sel <- inner_join(mdat_sel,details)
mdat_sel$BIRTH_DAY <- as.character(mdat_sel$BIRTH_DAY)

##################
##################

# check distribution of breeds and bdays across cohorts to choose a breed and bdays to form groups 

distribution_check <- mdat_sel %>%
  dplyr::select(breed,BIRTH_DAY,Cohort,collection_date)

# save distribution table 
distribution_check$grouping <- paste0(distribution_check$breed,"_",distribution_check$BIRTH_DAY)
x <- ctable(x = distribution_check$Cohort, y = distribution_check$grouping, prop = "r")
view(x, file = "distribution_check.html")

# print distribution table on console 
distribution_check <- xtabs(~ breed+BIRTH_DAY+Cohort, data=distribution_check)
ftable(distribution_check)

##################
##################

# decision : 
# group_A: breed: Duroc x Landrace, BIRTH_DAY: 2017-01-09 , 2017-01-10, 2017-01-11
# group_B: breed: Duroc x Large white, BIRTH_DAY: 2017-01-09 , 2017-01-10, 2017-01-11

##################
##################

# prepare metadata selections for groups A and B: 

group_A <- mdat_sel %>%
  filter(breed=="Duroc x Landrace") %>%
  filter(BIRTH_DAY=="2017-01-09"|BIRTH_DAY=="2017-01-10"|BIRTH_DAY=="2017-01-11")

group_B <- mdat_sel %>%
  filter(breed=="Duroc x Large white") %>%
  filter(BIRTH_DAY=="2017-01-09"|BIRTH_DAY=="2017-01-10"|BIRTH_DAY=="2017-01-11")


# dates selection for group A: 

group_A_Ja31 <- group_A %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

group_A_Fe7 <- group_A %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

group_A_Fe14 <- group_A %>% 
  filter(collection_date == "2017-02-14")

group_A_Fe21 <- group_A %>% 
  filter(collection_date == "2017-02-21")

group_A_Fe28 <- group_A %>% 
  filter(collection_date == "2017-02-28")

group_A_Ma3 <- group_A %>% 
  filter(collection_date == "2017-03-03")

A_a <- as.character(group_A_Ja31$ID)
A_b <- as.character(group_A_Fe7$ID)
A_c <- as.character(group_A_Fe14$ID)
A_d <- as.character(group_A_Fe21$ID)
A_e <- as.character(group_A_Fe28$ID)
A_f <- as.character(group_A_Ma3$ID)

writeLines(unlist(A_a), "guppy_input/piggies_group_A_Ja31_sel.txt", sep = " ")
writeLines(unlist(A_b), "guppy_input/piggies_group_A_Fe7_sel.txt", sep = " ")
writeLines(unlist(A_c), "guppy_input/piggies_group_A_Fe14_sel.txt", sep = " ")
writeLines(unlist(A_d), "guppy_input/piggies_group_A_Fe21_sel.txt", sep = " ")
writeLines(unlist(A_e), "guppy_input/piggies_group_A_Fe28_sel.txt", sep = " ")
writeLines(unlist(A_f), "guppy_input/piggies_group_A_Ma3_sel.txt", sep = " ")

# dates selection for group A: 

group_B_Ja31 <- group_B %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

group_B_Fe7 <- group_B %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

group_B_Fe14 <- group_B %>% 
  filter(collection_date == "2017-02-14")

group_B_Fe21 <- group_B %>% 
  filter(collection_date == "2017-02-21")

group_B_Fe28 <- group_B %>% 
  filter(collection_date == "2017-02-28")

group_B_Ma3 <- group_B %>% 
  filter(collection_date == "2017-03-03")


B_a <- as.character(group_B_Ja31$ID)
B_b <- as.character(group_B_Fe7$ID)
B_c <- as.character(group_B_Fe14$ID)
B_d <- as.character(group_B_Fe21$ID)
B_e <- as.character(group_B_Fe28$ID)
B_f <- as.character(group_B_Ma3$ID)


writeLines(unlist(B_a), "guppy_input/piggies_group_B_Ja31_sel.txt", sep = " ")
writeLines(unlist(B_b), "guppy_input/piggies_group_B_Fe7_sel.txt", sep = " ")
writeLines(unlist(B_c), "guppy_input/piggies_group_B_Fe14_sel.txt", sep = " ")
writeLines(unlist(B_d), "guppy_input/piggies_group_B_Fe21_sel.txt", sep = " ")
writeLines(unlist(B_e), "guppy_input/piggies_group_B_Fe28_sel.txt", sep = " ")
writeLines(unlist(B_f), "guppy_input/piggies_group_B_Ma3_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and NEO ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="Neomycin") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_Ja31 <- mdat_sel %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

mdat_Fe7 <- mdat_sel %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

mdat_Fe14 <- mdat_sel %>% 
  filter(collection_date == "2017-02-14")

mdat_Fe21 <- mdat_sel %>% 
  filter(collection_date == "2017-02-21")

mdat_Fe28 <- mdat_sel %>% 
  filter(collection_date == "2017-02-28")

mdat_Ma3 <- mdat_sel %>% 
  filter(collection_date == "2017-03-03")

a <- as.character(mdat_Ja31$ID)
b <- as.character(mdat_Fe7$ID)
c <- as.character(mdat_Fe14$ID)
d <- as.character(mdat_Fe21$ID)
e <- as.character(mdat_Fe28$ID)
f <- as.character(mdat_Ma3$ID)

writeLines(unlist(a), "guppy_input/piggies_CTRLNEO_Ja31_sel.txt", sep = " ")
writeLines(unlist(b), "guppy_input/piggies_CTRLNEO_Fe7_sel.txt", sep = " ")
writeLines(unlist(c), "guppy_input/piggies_CTRLNEO_Fe14_sel.txt", sep = " ")
writeLines(unlist(d), "guppy_input/piggies_CTRLNEO_Fe21_sel.txt", sep = " ")
writeLines(unlist(e), "guppy_input/piggies_CTRLNEO_Fe28_sel.txt", sep = " ")
writeLines(unlist(f), "guppy_input/piggies_CTRLNEO_Ma3_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : NEO and NEO+D ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Neomycin"|Cohort=="Neomycin+D-scour") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_Ja31 <- mdat_sel %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

mdat_Fe7 <- mdat_sel %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

mdat_Fe14 <- mdat_sel %>% 
  filter(collection_date == "2017-02-14")

mdat_Fe21 <- mdat_sel %>% 
  filter(collection_date == "2017-02-21")

mdat_Fe28 <- mdat_sel %>% 
  filter(collection_date == "2017-02-28")

mdat_Ma3 <- mdat_sel %>% 
  filter(collection_date == "2017-03-03")

a <- as.character(mdat_Ja31$ID)
b <- as.character(mdat_Fe7$ID)
c <- as.character(mdat_Fe14$ID)
d <- as.character(mdat_Fe21$ID)
e <- as.character(mdat_Fe28$ID)
f <- as.character(mdat_Ma3$ID)

writeLines(unlist(a), "guppy_input/piggies_NEONEOD_Ja31_sel.txt", sep = " ")
writeLines(unlist(b), "guppy_input/piggies_NEONEOD_Fe7_sel.txt", sep = " ")
writeLines(unlist(c), "guppy_input/piggies_NEONEOD_Fe14_sel.txt", sep = " ")
writeLines(unlist(d), "guppy_input/piggies_NEONEOD_Fe21_sel.txt", sep = " ")
writeLines(unlist(e), "guppy_input/piggies_NEONEOD_Fe28_sel.txt", sep = " ")
writeLines(unlist(f), "guppy_input/piggies_NEONEOD_Ma3_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : NEO and NEO+C ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Neomycin"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_Ja31 <- mdat_sel %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

mdat_Fe7 <- mdat_sel %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

mdat_Fe14 <- mdat_sel %>% 
  filter(collection_date == "2017-02-14")

mdat_Fe21 <- mdat_sel %>% 
  filter(collection_date == "2017-02-21")

mdat_Fe28 <- mdat_sel %>% 
  filter(collection_date == "2017-02-28")

mdat_Ma3 <- mdat_sel %>% 
  filter(collection_date == "2017-03-03")

a <- as.character(mdat_Ja31$ID)
b <- as.character(mdat_Fe7$ID)
c <- as.character(mdat_Fe14$ID)
d <- as.character(mdat_Fe21$ID)
e <- as.character(mdat_Fe28$ID)
f <- as.character(mdat_Ma3$ID)

writeLines(unlist(a), "guppy_input/piggies_NEONEOC_Ja31_sel.txt", sep = " ")
writeLines(unlist(b), "guppy_input/piggies_NEONEOC_Fe7_sel.txt", sep = " ")
writeLines(unlist(c), "guppy_input/piggies_NEONEOC_Fe14_sel.txt", sep = " ")
writeLines(unlist(d), "guppy_input/piggies_NEONEOC_Fe21_sel.txt", sep = " ")
writeLines(unlist(e), "guppy_input/piggies_NEONEOC_Fe28_sel.txt", sep = " ")
writeLines(unlist(f), "guppy_input/piggies_NEONEOC_Ma3_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and Ds ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_Ja31 <- mdat_sel %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

mdat_Fe7 <- mdat_sel %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

mdat_Fe14 <- mdat_sel %>% 
  filter(collection_date == "2017-02-14")

mdat_Fe21 <- mdat_sel %>% 
  filter(collection_date == "2017-02-21")

mdat_Fe28 <- mdat_sel %>% 
  filter(collection_date == "2017-02-28")

mdat_Ma3 <- mdat_sel %>% 
  filter(collection_date == "2017-03-03")

a <- as.character(mdat_Ja31$ID)
b <- as.character(mdat_Fe7$ID)
c <- as.character(mdat_Fe14$ID)
d <- as.character(mdat_Fe21$ID)
e <- as.character(mdat_Fe28$ID)
f <- as.character(mdat_Ma3$ID)

writeLines(unlist(a), "guppy_input/piggies_CTRLDs_Ja31_sel.txt", sep = " ")
writeLines(unlist(b), "guppy_input/piggies_CTRLDs_Fe7_sel.txt", sep = " ")
writeLines(unlist(c), "guppy_input/piggies_CTRLDs_Fe14_sel.txt", sep = " ")
writeLines(unlist(d), "guppy_input/piggies_CTRLDs_Fe21_sel.txt", sep = " ")
writeLines(unlist(e), "guppy_input/piggies_CTRLDs_Fe28_sel.txt", sep = " ")
writeLines(unlist(f), "guppy_input/piggies_CTRLDs_Ma3_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and Co ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_Ja31 <- mdat_sel %>% 
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01")

mdat_Fe7 <- mdat_sel %>% 
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07")

mdat_Fe14 <- mdat_sel %>% 
  filter(collection_date == "2017-02-14")

mdat_Fe21 <- mdat_sel %>% 
  filter(collection_date == "2017-02-21")

mdat_Fe28 <- mdat_sel %>% 
  filter(collection_date == "2017-02-28")

mdat_Ma3 <- mdat_sel %>% 
  filter(collection_date == "2017-03-03")

a <- as.character(mdat_Ja31$ID)
b <- as.character(mdat_Fe7$ID)
c <- as.character(mdat_Fe14$ID)
d <- as.character(mdat_Fe21$ID)
e <- as.character(mdat_Fe28$ID)
f <- as.character(mdat_Ma3$ID)

writeLines(unlist(a), "guppy_input/piggies_CTRLC_Ja31_sel.txt", sep = " ")
writeLines(unlist(b), "guppy_input/piggies_CTRLC_Fe7_sel.txt", sep = " ")
writeLines(unlist(c), "guppy_input/piggies_CTRLC_Fe14_sel.txt", sep = " ")
writeLines(unlist(d), "guppy_input/piggies_CTRLC_Fe21_sel.txt", sep = " ")
writeLines(unlist(e), "guppy_input/piggies_CTRLC_Fe28_sel.txt", sep = " ")
writeLines(unlist(f), "guppy_input/piggies_CTRLC_Ma3_sel.txt", sep = " ")


###########################################################################################
###########################################################################################


# these are going to guppy: 

# 1 piggies_sel.txt
# 1 pos_controls_sel.txt
# 6 piggies_*.txt
# 6 piggies_group_A_*.txt
# 6 piggies_group_B_*.txt
# 6 piggies_CTRLNEO_*.txt
# 6 piggies_NEONEOD_*.txt
# 6 piggies_NEONEOC_*.txt
# 6 piggies_CTRLDs_*.txt
# 6 piggies_CTRLC_*.txt

###########################################################################################
###########################################################################################


# NEXT: 

# move all the *_sel.txt files on to HPC here /shared/homes/12705859/phylosift_metapigs_20200225

# necessary input files in the sample directory: 
# all the .jplace files are 960 :
# (base) phylosift_metapigs_20200225/$ ls *.jplace | wc -l
# 960


# script to run: 

# nano guppy_epca.sh

# ```
#!/bin/bash
#PBS -l ncpus=4
#PBS -l walltime=48:00:00
#PBS -l mem=24g
# export OMP_NUM_THREADS=4
# cd $PBS_O_WORKDIR
# for f in /shared/homes/12705859/phylosift_metapigs_20200225/*.txt; do N=$(basename $f); /shared/homes/12705859/phylosift_v1.0.1/bin/guppy epca --prefix $N `cat $N`; done# ```
# ```
# works! 



# output files are: (extensions .xml, .edgediff, .trans, .jplace)
# move all the *_sel.txt.jplace and all the *_sel.txt.xml files to local machine /Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output
# .xml and .jplace files will be processed  with guppy_output_process.R 

###########################################################################################
###########################################################################################




