sample_df

sample_df$new_name <- paste0(sample_df$date,"_",sample_df$cohort,"_",sample_df$pig)
rownames(sample_df) <- sample_df[,5]

df
head(df)

df <- merge(no_reps_all, gtdbtk_bins_completeTaxa, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df)[colnames(df) == 'node'] <- 'gOTU'
df$gOTU <- paste0("gOTU_",df$gOTU)

NROW(unique(df$gOTU))
NROW(df)





# gOTU  

# columns to be kept 
df1 <- df
keep <- c("cohort","pig","bin","date","value","family")
df1 <- df1[ , (names(df1) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"

NROW(df1)

# assign a unique sample name 
df1$sample <- paste0(df1$date,"_",df1$cohort,"_",df1$pig)
# remove now pig and date (redundant)
df1$pig <- NULL
df1$date <- NULL
df1$cohort <- NULL
df1$bin <- NULL

# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(sample,family) %>%
  dplyr::summarise(all_bins_value = sum(value))

# long to wide 
df3 <- df2 %>%
  pivot_wider(names_from = sample, values_from = all_bins_value, values_fill = list(all_bins_value = 0)) 

# to matrix 
df3 <- as.data.frame(df3)
head(df3)

df3 <- df3[!is.na(df3[, 1]), ]

rownames(df3) <- df3$family

df3$family <- NULL





library(SIAMCAT)
#these are the counts and metadata from siamcat as they should look 
head(feat.crc.zeller)
head(meta.crc.zeller)
class(feat.crc.zeller)
class(meta.crc.zeller)

head(df3)
head(sample_df)




# turn into relative abundances - range 0 to 1 - (sample size normalization)
scale <- function(x, na.rm = FALSE) (x / sum(x))
normalized <- df3 %>%
  mutate_all(scale)




colSums(normalized)
head(normalized)
rownames(normalized) <- rownames(df3) 

# create metadata
empty_df = data.frame(
  sample = character())
df = rbind(
  empty_df,
  data.frame(
    sample = colnames(normalized)
  )
)

splitcols <- cSplit(df, "sample", "_")
aaa <- unite(splitcols,sampleID, 1:3, sep = "_", remove = FALSE)
names(aaa) <- c("sampleID", "date", "cohort", "pig")

df <- as.data.frame(aaa)
metadata <- data.frame(df[,-1], row.names=df[,1])
metadata$group <- paste0(metadata$date,"_",metadata$cohort)
label.normalized <- create.label(meta=metadata,
                                 label='group', 
                                 case='t4_NeoD',
                                 control='t4_Neomycin')
siamcat <- siamcat(feat=normalized,
                   label=label.normalized,
                   meta=metadata)
siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)

siamcat <- check.associations(
  siamcat,
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
model.evaluation.plot(siamcat)

model.interpretation.plot(
  siamcat,
  fn.plot = 'interpretation.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore')


####################################################################################################################################

