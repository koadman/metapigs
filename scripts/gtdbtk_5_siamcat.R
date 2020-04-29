



library(readr)
library(splitstackshape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(SIAMCAT)
library(matrixStats)
library(data.table)
library(pheatmap)


setwd("~/Desktop/metapigs_dry/gtdbtk")
basedir = "~/Desktop/metapigs_dry/"


# input files: 
# gtdbtk_bins_completeTaxa
# no_reps_all.csv (BINS COUNTS)
# dictionary

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

###### run this part if you want to only retain piglets that
# where present throughput trial (not euthanised)

# to_keep <- no_reps_all %>%
#   filter(date=="t8") %>%
#   dplyr::select(pig) %>%
#   distinct()
# 
# to_keep <- as.character(to_keep$pig)
# 
# no_reps_all <- subset(no_reps_all, (pig %in% to_keep))
# NROW(unique(no_reps_all$pig))

############################################################


######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv("gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################

# load gtdbtk dictionary

dict <- read_csv("bac120_arc122_dictionary",
                 col_types = cols(node = col_character()))

colnames(dict)[colnames(dict)=="node"] <- "gOTU"
dict$species <- gsub(" ","_", dict$species)

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

df$gOTU <- gsub("gOTU_","", df$gOTU)

NROW(df)
head(df)


df$gOTU <- paste0(df$phylum,"__",df$species,"__",df$gOTU)

df <- df %>% 
  dplyr::select(pig,bin,date,cohort,secondary_cluster,value,gOTU)

unique(df$gOTU)


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
  group_by(pig,date,gOTU) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df2)

# take the mean of each species by date:
df3 <- df2 %>%
  group_by(pig,gOTU,date) %>%
  dplyr::summarize(mean = mean(indiv_sum, na.rm=TRUE)) %>%
  mutate(perc=mean) %>%
  dplyr::select(pig,gOTU,date,perc) 
head(df2)


df3 <- as.data.frame(df3)
df3$sample = paste0(df3$date,"_",df3$pig)

# pivot wider
df3 <- df3 %>%
  dplyr::select(sample,gOTU,perc) %>%
  pivot_wider(names_from = sample, values_from = perc, values_fill = list(perc = 0))

feat <- as.data.frame(df3)
which(is.na(feat[,1]))


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

# create function to compare time points 

comparetimepoints_full <- function(t1,t2) { # where c is the cohort of interest, t0 is the start time point, t1 is the next time point)
  
  label.normalized <- create.label(meta=meta,
                                   label='date', 
                                   case= paste0(t2),
                                   control= paste0(t1))
  
  siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
  siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
  
  # check for significant associations
  siamcat <- check.associations(
    siamcat,
    fn.plot = paste0("gt_siamcatA_",t1,t2,".pdf"),
    sort.by = 'fc',
    alpha = 0.05, 
    mult.corr = "fdr",
    detect.lim = 10 ^-6,
    plot.type = "quantile.box",
    panels = c("fc", "prevalence", "auroc"))
  
  
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
    fn.plot = paste0("gt_siamcatH_",t1,t2,".pdf"),
    consens.thres = 0.5, 
    limits = c(-3, 3),
    heatmap.type = 'zscore')
  
}


# run the timepoint comparisons you want: (outputs are two plots each, automatically saved)
# 1 week interval: 
comparetimepoints_full("t0","t2")
comparetimepoints_full("t2","t4")
comparetimepoints_full("t4","t6")
comparetimepoints_full("t6","t8")
comparetimepoints_full("t8","t10")
# 2 weeks interval:
comparetimepoints_full("t0","t4")
comparetimepoints_full("t4","t8")
# 4 weeks interval:
comparetimepoints_full("t0","t8")


############################################################################################################################################
############################################################################################################################################


# comparing cohorts 

######################################################################
######################################################################

# create function to compare cohorts: 

# prepare empty df to be filled
dfz <- data.frame(
  fc = numeric(),
  p.adj = numeric(),
  gOTU = character(),
  comparison = character(),
  stringsAsFactors = FALSE
)

# function to compare cohorts: 
comparebetween <- function(t,c2,c1) { 
  
  label.normalized <- create.label(meta=meta,
                                   label='group', 
                                   case= paste0(t,"_",c2),
                                   control= paste0(t,"_",c1))
  
  siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
  siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
  
  # check for significant associations
  siamcat <- check.associations(
    siamcat,
    sort.by = 'fc',
    alpha = 0.05,
    mult.corr = "fdr",
    detect.lim = 10 ^-6,
    prompt = FALSE,
    panels = c("fc", "prevalence", "auroc"))
  
  yy <- SIAMCAT::associations(siamcat)
  yy$gOTU <- rownames(yy)
  yy <- yy %>% filter(p.adj<0.05)
  
  
  if (NROW(yy)>0) {
    yy$comparison=paste0(t,"_",c2,"_",c1)
    yy$feat <- rownames(yy)
    yy <- yy %>%
      dplyr::select(comparison,gOTU,fc,p.adj)
    dfz <- rbind(dfz,yy)
  } else {
    print ("No significant associations found. No plot will be produced.")
  }

  return(dfz)
}


NROW(last)

# run the association tests and catch all the df outputs into one single df
last <- rbind(comparebetween("t0","DScour","Control"),
              comparebetween("t2","DScour","Control"),
              comparebetween("t4","DScour","Control"),
              comparebetween("t6","DScour","Control"),
              comparebetween("t8","DScour","Control"),
              comparebetween("t10","DScour","Control"),
              comparebetween("t0","ColiGuard","Control"),
              comparebetween("t2","ColiGuard","Control"),
              comparebetween("t4","ColiGuard","Control"),
              comparebetween("t6","ColiGuard","Control"),
              comparebetween("t8","ColiGuard","Control"),
              comparebetween("t10","ColiGuard","Control"),
              comparebetween("t0","Neomycin","Control"),
              comparebetween("t2","Neomycin","Control"),
              comparebetween("t4","Neomycin","Control"),
              comparebetween("t6","Neomycin","Control"),
              comparebetween("t8","Neomycin","Control"),
              comparebetween("t10","Neomycin","Control"),
              comparebetween("t0","NeoD","Neomycin"),
              comparebetween("t2","NeoD","Neomycin"),
              comparebetween("t4","NeoD","Neomycin"),
              comparebetween("t6","NeoD","Neomycin"),
              comparebetween("t8","NeoD","Neomycin"),
              comparebetween("t10","NeoD","Neomycin"),
              comparebetween("t0","NeoC","Neomycin"),
              comparebetween("t2","NeoC","Neomycin"),
              comparebetween("t4","NeoC","Neomycin"),
              comparebetween("t6","NeoC","Neomycin"),
              comparebetween("t8","NeoC","Neomycin"),
              comparebetween("t10","NeoC","Neomycin"))


last$p.adj <- NULL


sink(file = "gt_siamcat_cohortsbetween.txt", 
     append = FALSE, type = c("output"))
last
sink()


######################################################################
######################################################################

# Function to plot the significant changes per cohort for each time interval : 

comparewithin_plots <- function(c,t1,t2) { # where c is the cohort of interest, t0 is the start time point, t1 is the next time point)
  
  label.normalized <- create.label(meta=meta,
                                   label='group', 
                                   case= paste0(t2,"_",c),
                                   control= paste0(t1,"_",c))
  
  siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
  siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
  
  
  # check for significant associations
  siamcat <- check.associations(
    siamcat,
    fn.plot = paste0("gt_siamcatA_",c,"_",t1,"_",t2,".pdf"),
    sort.by = 'fc',
    alpha = 0.05, 
    mult.corr = "fdr",
    detect.lim = 10 ^-6,
    plot.type = "quantile.box",
    panels = c("fc", "prevalence", "auroc"))
  
  
}


# Plot : 

comparewithin_plots("Control","t0","t2")
comparewithin_plots("Control","t2","t4")
comparewithin_plots("Control","t4","t6")
comparewithin_plots("Control","t6","t8")
comparewithin_plots("Control","t8","t10")
comparewithin_plots("Control","t4","t8")
comparewithin_plots("Control","t4","t10")
comparewithin_plots("Control","t2","t8")
#
comparewithin_plots("DScour","t0","t2")
comparewithin_plots("DScour","t2","t4")
comparewithin_plots("DScour","t4","t6")
comparewithin_plots("DScour","t6","t8")
comparewithin_plots("DScour","t8","t10")
comparewithin_plots("DScour","t4","t8")
comparewithin_plots("DScour","t4","t10")
comparewithin_plots("DScour","t2","t8")
#
comparewithin_plots("ColiGuard","t0","t2")
comparewithin_plots("ColiGuard","t2","t4")
comparewithin_plots("ColiGuard","t4","t6")
comparewithin_plots("ColiGuard","t6","t8")
comparewithin_plots("ColiGuard","t8","t10")
comparewithin_plots("ColiGuard","t4","t8")
comparewithin_plots("ColiGuard","t4","t10")
comparewithin_plots("ColiGuard","t2","t8")
#
comparewithin_plots("Neomycin","t0","t2")
comparewithin_plots("Neomycin","t2","t4")
comparewithin_plots("Neomycin","t4","t6")
comparewithin_plots("Neomycin","t6","t8")
comparewithin_plots("Neomycin","t8","t10")
comparewithin_plots("Neomycin","t4","t8")
comparewithin_plots("Neomycin","t4","t10")
comparewithin_plots("Neomycin","t2","t8")
#
comparewithin_plots("NeoD","t0","t2")
comparewithin_plots("NeoD","t2","t4")
comparewithin_plots("NeoD","t4","t6")
comparewithin_plots("NeoD","t6","t8")
comparewithin_plots("NeoD","t8","t10")
comparewithin_plots("NeoD","t4","t8")
comparewithin_plots("NeoD","t4","t10")
comparewithin_plots("NeoD","t2","t8")
#
comparewithin_plots("NeoC","t0","t2")
comparewithin_plots("NeoC","t2","t4")
comparewithin_plots("NeoC","t4","t6")
comparewithin_plots("NeoC","t6","t8")
comparewithin_plots("NeoC","t8","t10")
comparewithin_plots("NeoC","t4","t8")
comparewithin_plots("NeoC","t4","t10")
comparewithin_plots("NeoC","t2","t8")


######################################################################

# Function to obtain the (data only) significant changes per cohort for each time interval : 

dfz <- data.frame(
  fc = numeric(),
  p.adj = numeric(),
  gOTU = character(),
  comparison = character(),
  stringsAsFactors = FALSE
)

comparewithin <- function(c,t1,t2) { # where c is the cohort of interest, t0 is the start time point, t1 is the next time point)
  
  label.normalized <- create.label(meta=meta,
                                   label='group', 
                                   case= paste0(t2,"_",c),
                                   control= paste0(t1,"_",c))
  
  siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
  siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
  
  # check for significant associations
  siamcat <- check.associations(
    siamcat,
    sort.by = 'fc',
    alpha = 0.05,
    mult.corr = "fdr",
    detect.lim = 10 ^-6,
    prompt = FALSE,
    panels = c("fc", "prevalence", "auroc"))
  
  yy <- SIAMCAT::associations(siamcat)
  yy$gOTU <- rownames(yy)
  yy <- yy %>% filter(p.adj<0.05)
  
  
  if (NROW(yy)>0) {
    yy$comparison=paste0(c,"_",t1,"_",t2)
    yy$feat <- rownames(yy)
    yy <- yy %>%
      dplyr::select(comparison,gOTU,fc,p.adj)
    dfz <- rbind(dfz,yy)
  } else {
    print ("No significant associations found. No plot will be produced.")
  }
  
  return(dfz)
}



last <- rbind(comparewithin("Control","t0","t2"),
              comparewithin("Control","t2","t4"),
              comparewithin("Control","t4","t6"),
              comparewithin("Control","t6","t8"),
              comparewithin("Control","t8","t10"),
              comparewithin("Control","t4","t8"),
              comparewithin("Control","t4","t10"),
              comparewithin("Control","t2","t8"),
              #
              comparewithin("DScour","t0","t2"),
              comparewithin("DScour","t2","t4"),
              comparewithin("DScour","t4","t6"),
              comparewithin("DScour","t6","t8"),
              comparewithin("DScour","t8","t10"),
              comparewithin("DScour","t4","t8"),
              comparewithin("DScour","t4","t10"),
              comparewithin("DScour","t2","t8"),
              #
              comparewithin("ColiGuard","t0","t2"),
              comparewithin("ColiGuard","t2","t4"),
              comparewithin("ColiGuard","t4","t6"),
              comparewithin("ColiGuard","t6","t8"),
              comparewithin("ColiGuard","t8","t10"),
              comparewithin("ColiGuard","t4","t8"),
              comparewithin("ColiGuard","t4","t10"),
              comparewithin("ColiGuard","t2","t8"),
              #
              comparewithin("Neomycin","t0","t2"),
              comparewithin("Neomycin","t2","t4"),
              comparewithin("Neomycin","t4","t6"),
              comparewithin("Neomycin","t6","t8"),
              comparewithin("Neomycin","t8","t10"),
              comparewithin("Neomycin","t4","t8"),
              comparewithin("Neomycin","t4","t10"),
              comparewithin("Neomycin","t2","t8"),
              #
              comparewithin("NeoD","t0","t2"),
              comparewithin("NeoD","t2","t4"),
              comparewithin("NeoD","t4","t6"),
              comparewithin("NeoD","t6","t8"),
              comparewithin("NeoD","t8","t10"),
              comparewithin("NeoD","t4","t8"),
              comparewithin("NeoD","t4","t10"),
              comparewithin("NeoD","t2","t8"),
              #
              comparewithin("NeoC","t0","t2"),
              comparewithin("NeoC","t2","t4"),
              comparewithin("NeoC","t4","t6"),
              comparewithin("NeoC","t6","t8"),
              comparewithin("NeoC","t8","t10"),
              comparewithin("NeoC","t4","t8"),
              comparewithin("NeoC","t4","t10"),
              comparewithin("NeoC","t2","t8"))

last_to_use <- last 

NROW(last_to_use)
last_to_use$p.adj <- NULL

toplot <- last_to_use %>%
  dplyr::select(gOTU,comparison,fc)



# save hits as text file 
sink(file = "gt_siamcatA_cohortswithin.txt", 
     append = TRUE, type = c("output"))
last_to_use
sink()



toplot <- cSplit(toplot, "comparison", "_")
toplot$interval <- paste0(toplot$comparison_2,"_",toplot$comparison_3)

toplot <- toplot %>%
  dplyr::select(gOTU,comparison_1,interval,fc)

colnames(toplot) <- c("gOTU","cohort","interval","fc")


# I am adding this here because the pdf below often doesn t get printed
closeAllConnections()


# split df by time interval 
multiple_DFs <- split( toplot , f = toplot$interval )


pdf("gt_siamcatA_cohortswithin.pdf")
for (single_DF in multiple_DFs) {
  
  
  DF <- as.data.frame(single_DF)
  
  interval <- DF$interval[1]
  
  DF_wide <- DF %>%
    pivot_wider(names_from = cohort, values_from = fc, values_fill = list(fc = 0))
  
  if (NCOL(DF_wide) > 1 & NROW(DF_wide) >1 ) {
    
    DF_wide <- as.data.frame(DF_wide)
    
    rownames(DF_wide) <- DF_wide[,1]
    DF_wide[,1] <- NULL
    DF_wide$interval <- NULL
    
    m <- as.matrix(DF_wide)
    
    pheatmap(m, fontsize_row = 5, main = interval, display_numbers = TRUE, cluster_cols = FALSE)
    
  }
  
  
}
dev.off()





# matrix subsetting when heatmap too crowded : 

  
  # mat <- m
  # 
  # # relevant metrics per row
  # row_med <- matrixStats::rowMedians(mat)
  # row_vars <- matrixStats::rowVars(mat)
  # row_maxs <- matrixStats::rowMaxs(mat)
  # row_qntl90 <- matrixStats::rowQuantiles(mat, probs = 0.9)
  # 
  # # top 50% utility function
  # top5 <- function(x) {
  #   x >= quantile(x, 0.50)
  # }
  # 
  # # combine all conditions
  # row_idx <- top5(row_vars) & top5(row_maxs - row_med) & top5(row_qntl90 - row_med)
  # # subscript
  # n <- mat[row_idx, , drop = FALSE]

