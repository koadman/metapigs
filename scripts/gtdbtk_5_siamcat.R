



library(readr)
library(splitstackshape)
library(tidyr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(SIAMCAT)


setwd("~/Desktop/metapigs_dry/gtdbtk")
basedir = "~/Desktop/metapigs_dry/"


# input files: 
# gtdbtk_bins_completeTaxa
# no_reps_all.csv (BINS COUNTS)

# OUTPUTS:

# 


######################################################################

# a look at SIAMCAT dummy dataset:

#these are the counts and metadata from siamcat as they should look 
head(feat.crc.zeller)
head(meta.crc.zeller)
class(feat.crc.zeller)
class(meta.crc.zeller)

# the sum of each column (1 column: 1 sample) always equals 1. 
# this implies SIAMCAT needs ratios
sum(feat.crc.zeller[,2])

######################################################################

# counts data 

no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)


######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv("gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################

# merge info 

NROW(gtdbtk_bins)
NROW(no_reps_all)
head(gtdbtk_bins)
head(no_reps_all)
df <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df)[colnames(df) == 'node'] <- 'gOTU'
df$gOTU <- paste0("gOTU_",df$gOTU)

NROW(unique(df$gOTU))
NROW(df)

head(df)

# if we want to use below's code with gOTUs, there s no issues. 
# but if you wnat to select either family, species, genus, 
# you d run into problems of non-uniqueness. 
# way around it: add the gOTU_ID to each of the taxonomic levels 
df <- df %>%
  dplyr::mutate(species=paste0(species," ",gOTU)) %>%
  dplyr::mutate(genus=paste0(genus," ",gOTU)) %>%
  dplyr::mutate(family=paste0(family," ",gOTU)) %>%
  dplyr::mutate(order=paste0(order," ",gOTU)) %>%
  dplyr::mutate(class=paste0(class," ",gOTU)) %>%
  dplyr::mutate(phylum=paste0(phylum," ",gOTU)) %>%
  dplyr::mutate(domain=paste0(domain," ",gOTU))
  



######################################################################

# CREATE COUNTS TABLE (like feat.crc.zeller)

df1 <- df

# lib size normalization
df2 <- df1 %>% 
  group_by(pig,date) %>% 
  mutate(norm_value = (value/sum(value))) %>% 
  dplyr::select(-value)

head(df2)
colnames(df2)

# sum all the norm values that fall within same pig,date,species:
df2 <- df2 %>%  
  group_by(pig,date,species) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df2)

# take the mean of each species by date:
df3 <- df2 %>%
  group_by(pig,species,date) %>%
  dplyr::summarize(mean = mean(indiv_sum, na.rm=TRUE)) %>%
  mutate(perc=mean) %>%
  dplyr::select(pig,species,date,perc) 
head(df2)


df3 <- as.data.frame(df3)
df3$sample = paste0(df3$date,"_",df3$pig)

# pivot wider
df3 <- df3 %>%
  dplyr::select(sample,species,perc) %>%
  pivot_wider(names_from = sample, values_from = perc, values_fill = list(perc = 0))

feat <- as.data.frame(df3)

rownames(feat) <- feat[,1]
feat[,1] <- NULL


head(feat)
dim(feat)

# is the sum of each columns 1? 
sum(feat[,4])
# yes 

# ready! 

######################################################################

# CREATE METADATA TABLE (like meta.crc.zeller)

theseAREtheSamples <- as.data.frame(colnames(feat))
colnames(theseAREtheSamples) <- "sample"

df1$sample <- paste0(df1$date,"_",df1$pig)
df1 <- df1 %>%
  dplyr::select(sample,cohort,pig,date) %>%
  distinct()
# add anther grouping: date+cohort:
df1$group <- paste0(df1$date,"_",df1$cohort)

head(theseAREtheSamples)

meta <- left_join(theseAREtheSamples,df1)

rownames(meta) <- meta[,1]
meta[,1] <- NULL
# ready! 

######################################################################

# SIAMCAT starts! 

class(feat)
class(feat.crc.zeller)
class(meta)
class(meta.crc.zeller)

head(feat)
head(feat.crc.zeller)
head(meta)
head(meta.crc.zeller)

colnames(feat)==rownames(meta)





######################################################################
######################################################################
label.normalized <- create.label(meta=meta,
                                 label='date', 
                                 case='t2',
                                 control='t0')

siamcat <- siamcat(feat=feat,
                   label=label.normalized,
                   meta=meta)

siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  fn.plot = 'gt_siamcatA_t0t2.pdf',
  sort.by = 'fc',
  alpha = 0.05,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  panels = c("fc", "prevalence", "auroc"))

############################################

# Model building

siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)

siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 2
)

siamcat <- train.model(
  siamcat,
  method = "lasso"
)

siamcat <- make.predictions(siamcat)
pred_matrix <- pred_matrix(siamcat)
siamcat <-  evaluate.predictions(siamcat)
#model.evaluation.plot(siamcat)

model.interpretation.plot(
  siamcat,
  fn.plot = 'gt_siamcatH_t0t2.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore')

######################################################################
######################################################################

label.normalized <- create.label(meta=meta,
                                 label='date', 
                                 case='t4',
                                 control='t2')

siamcat <- siamcat(feat=feat,
                   label=label.normalized,
                   meta=meta)

siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  fn.plot = 'gt_siamcatA_t2t4.pdf',
  sort.by = 'fc',
  alpha = 0.05,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  panels = c("fc", "prevalence", "auroc"))

############################################

# Model building

siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)

siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 2
)

siamcat <- train.model(
  siamcat,
  method = "lasso"
)

siamcat <- make.predictions(siamcat)
pred_matrix <- pred_matrix(siamcat)
siamcat <-  evaluate.predictions(siamcat)
#model.evaluation.plot(siamcat)

model.interpretation.plot(
  siamcat,
  fn.plot = 'gt_siamcatH_t2t4.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore')

######################################################################
######################################################################


label.normalized <- create.label(meta=meta,
                                 label='date', 
                                 case='t6',
                                 control='t4')

siamcat <- siamcat(feat=feat,
                   label=label.normalized,
                   meta=meta)

siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  fn.plot = 'gt_siamcatA_t4t6.pdf',
  sort.by = 'fc',
  alpha = 0.05,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  panels = c("fc", "prevalence", "auroc"))

############################################

# Model building

siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)

siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 2
)

siamcat <- train.model(
  siamcat,
  method = "lasso"
)

siamcat <- make.predictions(siamcat)
pred_matrix <- pred_matrix(siamcat)
siamcat <-  evaluate.predictions(siamcat)
#model.evaluation.plot(siamcat)

model.interpretation.plot(
  siamcat,
  fn.plot = 'gt_siamcatH_t4t6.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore')

######################################################################
######################################################################


# una prova di trattamento! 


label.normalized <- create.label(meta=meta,
                                 label='group', 
                                 case='t4_NeoD',
                                 control='t4_Neomycin')

siamcat <- siamcat(feat=feat,
                   label=label.normalized,
                   meta=meta)

siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  fn.plot = 'gt_siamcatA_t4_NeoNeoD.pdf',
  sort.by = 'fc',
  alpha = 0.05,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  panels = c("fc", "prevalence", "auroc"))

############################################

# Model building

siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)

siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 2
)

siamcat <- train.model(
  siamcat,
  method = "lasso"
)

siamcat <- make.predictions(siamcat)
pred_matrix <- pred_matrix(siamcat)
siamcat <-  evaluate.predictions(siamcat)
#model.evaluation.plot(siamcat)

model.interpretation.plot(
  siamcat,
  fn.plot = 'gt_siamcatH_t4_NeoNeoD.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore')

######################################################################
######################################################################



# una prova di trattamento! 


label.normalized <- create.label(meta=meta,
                                 label='group', 
                                 case='t4_DScour',
                                 control='t4_Control')

siamcat <- siamcat(feat=feat,
                   label=label.normalized,
                   meta=meta)

siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  fn.plot = 'gt_siamcatA_t4_CtrlDScour.pdf',
  sort.by = 'fc',
  alpha = 0.05,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  panels = c("fc", "prevalence", "auroc"))

############################################

# Model building

siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)

siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 2
)

siamcat <- train.model(
  siamcat,
  method = "lasso"
)

siamcat <- make.predictions(siamcat)
pred_matrix <- pred_matrix(siamcat)
siamcat <-  evaluate.predictions(siamcat)
#model.evaluation.plot(siamcat)

model.interpretation.plot(
  siamcat,
  fn.plot = 'gt_siamcatH_t4_CtrlDScour.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore')

######################################################################
######################################################################





####################################################################################################################################

