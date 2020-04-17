
library(readr)
library(dplyr)
library(readxl)
library(splitstackshape)
library(data.table)
library(robCompositions)
library(tidyr)
library(ggbiplot)
library(magrittr)
library(ggpubr)



setwd("/Users/12705859/Desktop/metapigs_base/sortmerna")
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"



######################################################################

# upload sortmerna output (from the original output we retained columns:plate_well, 16Sgene ID, e-value)
# and we filtered out all the 16S rRNA genes below e-30 threshold 
so <- read_table2("sortmeall_evaluefiltered.tsv", col_names = FALSE)
so$X1 <- gsub("_S","", so$X1)
so <- so[,1:2]

so_work <- so

######################################################################

# load metadata 
mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)
colnames(mdat)

colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'

mdat <- mdat %>%
  dplyr::select(DNA_plate,DNA_well,Cohort,collection_date,isolation_source)

######################################################################


# further reducing size: AIM: 

so_work$count <- as.numeric(paste0(1))

NROW(so_work)
so_done <- setDT(so_work)[,.(A = sum(count)), by = 'X1,X2']
NROW(so_done)
# much an improvement! 

#############

# once size of sortme file is reduced, we can parse it: 
so_done_temp <- cSplit(so_done,"X1","_")
so_done_temp$DNA_plate <- paste0(so_done_temp$X1_1,"_",so_done_temp$X1_2)
so_done_temp$DNA_well <- paste0(so_done_temp$X1_3)
head(so_done_temp)

so_done_temp <- so_done_temp %>%
  dplyr::select(DNA_plate,DNA_well,X2,A)

colnames(so_done_temp)[colnames(so_done_temp) == 'X2'] <- 'rRNA16S'
colnames(so_done_temp)[colnames(so_done_temp) == 'A'] <- 'count'
head(so_done_temp)

###############################

# time to merge to metadata!

NROW(so_done_temp)
NROW(mdat)

df <- left_join(so_done_temp,mdat)
NROW(df)

head(df)


# tM <- "2017-01-30"
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-01-30", 
  replacement = "tM", 
  fixed = TRUE)
# t0 <- "2017-01-31" "2017-02-01" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-01-31", 
  replacement = "t0", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-01", 
  replacement = "t0", 
  fixed = TRUE)

# t1 <- "2017-02-03" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-03", 
  replacement = "t1", 
  fixed = TRUE)

# t2 <- "2017-02-06" "2017-02-07" "2017-02-08"
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-06", 
  replacement = "t2", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-07", 
  replacement = "t2", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-08", 
  replacement = "t2", 
  fixed = TRUE)

# t3 <- "2017-02-10" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-10", 
  replacement = "t3", 
  fixed = TRUE)

# t4 <- "2017-02-14"
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-14", 
  replacement = "t4", 
  fixed = TRUE)

# t4 <- "2017-02-16" "2017-02-17" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-16", 
  replacement = "t5", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-17", 
  replacement = "t5", 
  fixed = TRUE)

# t6 <- "2017-02-21" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-21", 
  replacement = "t6", 
  fixed = TRUE)

# t7 <- "2017-02-24" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-24", 
  replacement = "t7", 
  fixed = TRUE)

# t8 <- "2017-02-28" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-02-28", 
  replacement = "t8", 
  fixed = TRUE)

# t9 <- "2017-03-03" 
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-03-03", 
  replacement = "t9", 
  fixed = TRUE)

# t10 <- "2017-03-06" "2017-03-07" "2017-03-08" "2017-03-09" "2017-03-10"
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-03-06", 
  replacement = "t10", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-03-07", 
  replacement = "t10", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-03-08", 
  replacement = "t10", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-03-09", 
  replacement = "t10", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-03-10", 
  replacement = "t10", 
  fixed = TRUE)

df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2017-08-14", 
  replacement = "no-t-pos", 
  fixed = TRUE)
df[6] <- lapply(
  df[6], 
  gsub, 
  pattern = "2018-01-24", 
  replacement = "no-t-pos", 
  fixed = TRUE)

# no-t-neg for negative control
df <- df %>%
  mutate(collection_date = if_else(is.na(collection_date), "no-t-neg", collection_date))

unique(df$collection_date)

df$sample <- paste0(df$collection_date,"_",df$isolation_source)






# normalization for library size 
df1 <- df %>%
  dplyr::select(sample,rRNA16S,count) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(norm_count=count/sum(count)) 
NROW(df1)
head(df1)

# sum all the norm values that fall within same pig,date,taxa_2
df2 <- df1 %>%
  dplyr::group_by(sample,rRNA16S) %>%
  dplyr::summarise(last_count = sum(norm_count))
head(df2)


# filter out pos and neg controls (they add unwanted variance and we want to look at just piglet samples)
df2.5 <- df2 %>%
  filter(., !grepl("no-t",sample))  %>%
  filter(., !grepl("tM",sample)) 

# long to wide format
df3 <- df2.5 %>%
  pivot_wider(names_from = rRNA16S, values_from = last_count, values_fill = list(last_count = 0)) 
head(df3)

df_wide <- as.data.frame(df3)


dim(df_wide)
rownames(df_wide) <- df_wide[,1]
df_wide[,1] <- NULL

# it's 1, right! 
rowSums(df_wide)


x <- df_wide[,order(colSums(df_wide),decreasing=T)]




# using all rRNA16S genes 
mtcars.pca1 <- prcomp(x, center = TRUE,scale. = TRUE)
dates <- substr(rownames(df_wide), start = 1, stop = 3)  %<>%
  gsub('_$', '', .) # removes the last _
plot1 <- ggbiplot(mtcars.pca1, 
                 groups=dates, 
                 labels = NULL, 
                 var.axes=FALSE) +
  geom_point(aes(colour=dates), size = 0.01) +
  theme(panel.grid.major = element_blank(), 
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(), 
             axis.line = element_line(colour = "black"),
             legend.title = element_text()) +
  labs(colour="collection date") +
  ggtitle("All") + guides(colour = guide_legend(nrow = 1))
plot1
###


# using 150 most abundant rRNA16S genes 
mtcars.pca2 <- prcomp(x[,1:150], center = TRUE,scale. = TRUE)
dates <- substr(rownames(df_wide), start = 1, stop = 3)  %<>%
  gsub('_$', '', .) # removes the last _
plot2 <- ggbiplot(mtcars.pca2, 
         groups=dates, 
         labels = NULL,
         var.axes=FALSE) +
  geom_point(aes(colour=dates), size = 0.01) +
  theme(panel.grid.major = element_blank(), 
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(), 
             axis.line = element_line(colour = "black"),
             legend.title = element_text()) +
  labs(colour="collection date") +
  ggtitle("150 most abundant") + guides(colour = guide_legend(nrow = 1))
plot2
###



# using 20 most abundant rRNA16S genes 
mtcars.pca3 <- prcomp(x[,1:20], center = TRUE,scale. = TRUE)
dates <- substr(rownames(df_wide), start = 1, stop = 3)  %<>%
  gsub('_$', '', .) # removes the last _
plot3 <- ggbiplot(mtcars.pca3, 
                 groups=dates, 
                 labels = NULL,
                 var.axes=FALSE) +
  geom_point(aes(colour=dates), size = 0.01) +
  theme(panel.grid.major = element_blank(), 
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(), 
             axis.line = element_line(colour = "black"),
             legend.title = element_text()) +
  labs(colour="collection date")  +
  ggtitle("20 most abundant") + guides(colour = guide_legend(nrow = 1))
plot3
###



tosave <- ggarrange(plot1, plot2, plot3,
                    ncol = 3, 
                    labels=c("A","B","C"), 
                    common.legend = TRUE, 
                    align = "h", legend = "top",
                    widths = c(2,2,2),
                    heights = c(2,2,2))

pdf("sortmerna_all.pdf", width = 8, height = 5)
tosave
dev.off()




