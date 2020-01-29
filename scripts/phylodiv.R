
#  R version 3.6.1 (2019-07-05)

# 0   # loading input data
# 1   # batch effect (plot,removal and plot again)
# 1.1 # merge alpha div data with metadata -> boggo
# 1.2 # merge beta div data with metadata -> coggo
# 2   # plot samples distribution (plate vs time; cohorts vs time)
# 2.1 # boggo and coggo: formatting & averaging duplicates
# 3   # DELTAS
# 4   # plot ALPHA diversity (all timepoints)
# 5   # plot BETA diversity (all timepoints)
# 6   # plot ALPHA diversity (at pig trial start) 
# 7   # plot BETA diversity (at pig trial start) 
# 8   # plot distribution of breeds and bdays among cohorts + cohorts distr plates 
# 9   # p-values
# 10 # plot p-values (starting factors)
# 11  # prepare input files for guppy (select by time interval) and plot output 

###########################################################################################

# install packages

#pkgs <- c("ggbiplot","ggpubr","sva","tidyverse","broom","cowplot","data.table","dunn.test","plyr",
#           "dplyr","forcats","ggplot2","gridExtra","plotrix","readr","readxl","tidyr","varhandle","tibble","purr","remotes")

#install.packages(pkgs[], repos='https://cran.rstudio.com', dependencies = TRUE)  
#remotes::install_github("vqv/ggbiplot")

###########################################################################################

#BiocManager::install("sva")
#BiocManager::install("genefilter")

# load libraries
library(magick)
library(tiff)
library(rstatix)
library(ggbiplot) # ggbiplot_0.55 
library(ggpubr) # ggpubr_0.2.4 
library(sva) # sva_3.32.1  
library(tidyverse) # tidyverse_1.3.0 
library(broom) # broom_0.5.2
library(cowplot) # cowplot_1.0.0
library(data.table) # data.table_1.12.8  
library(dunn.test) # dunn.test_1.3.5 
library(plyr) # plyr_1.8.5 
library(dplyr) # dplyr_0.8.3 
library(forcats) # forcats_0.4.0 
library(ggplot2) # ggplot2_3.2.1
library(gridExtra) # gridExtra_2.3     
library(plotrix) # plotrix_3.7-7
library(readr) # readr_1.3.1
library(readxl) # readxl_1.3.1
library(tidyr) # tidyr_1.0.0
library(varhandle) # varhandle_2.0.4
library(tibble) # tibble_2.1.3 
library(purrr) # purrr_0.3.3
library(openxlsx)
library(genefilter)
library(compareGroups)
library(splitstackshape)

# from local 
setwd("/Users/12705859/Desktop/metapigs_base/phylosift")
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"

# from HPC
#setwd("/shared/homes/s1/pig_microbiome/metapigs_base")
#basedir = "/shared/homes/s1/pig_microbiome/metapigs_base/input_files/"

###########################################################################################

# 0   # loading input data

# tiffs (timelines)
timeline <- image_read(paste0(basedir,"Slide01.tiff"))
timeline_31_Jan <- image_read(paste0(basedir,"Slide02.tiff"))
timeline_3_Feb <- image_read(paste0(basedir,"Slide03.tiff"))
timeline_7_Feb <- image_read(paste0(basedir,"Slide04.tiff"))
timeline_10_Feb <- image_read(paste0(basedir,"Slide05.tiff"))
timeline_14_Feb <- image_read(paste0(basedir,"Slide06.tiff"))
timeline_17_21_Feb <- image_read(paste0(basedir,"Slide07.tiff"))
timeline_24_28_Feb <- image_read(paste0(basedir,"Slide08.tiff"))
timeline_3_10_Mar <- image_read(paste0(basedir,"Slide09.tiff"))
timeline_deltas <- image_read(paste0(basedir,"Slide10.tiff"))
timeline_deltas_unroo <- image_read(paste0(basedir,"Slide11.tiff"))
timeline_deltas_bw <- image_read(paste0(basedir,"Slide12.tiff"))

# load metadata 
mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

# load details (breed, line, bday, mothers)
details <- read_excel(paste0(basedir, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'pig'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'stig'
colnames(details)[colnames(details) == '...8'] <- 'breed'
details$pig <- gsub("G","",details$pig)
details$pig <- gsub("T","",details$pig)

details <- details %>%
  select(pig,BIRTH_DAY,LINE,breed,stig,nurse)

# formatting metadata column names 
mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == '*sample_name'] <- 'sample_name'

# load alpha pd
fpddat<-read.table(paste0(basedir,"fpdalpha_div.tsv"),header=T,stringsAsFactors=F)
# load beta div
pcadat<-read.table(paste0(basedir,"new.proj"),header=T,stringsAsFactors=F)

###########################################################################################

# create workbook to add stats 

wb <- createWorkbook()

###########################################################################################


# 1   # batch effect

# Plots the batch effect (both alpha and beta div)

# as pcadat has duplicate rows, keep unique rows: 
pcadat <- unique(pcadat)

fpddat$sid <- NULL
pcadat$s_id <- NULL

# re-order plates
fpddat$DNA_plate <- factor(fpddat$DNA_plate, 
                           levels=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))
pcadat$DNA_plate <- factor(pcadat$DNA_plate, 
                           levels=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))

palette <- c("black","red","green3","blue","cyan","magenta","yellow","gray","orange","brown")

pdf("out/batch_pca.pdf")
par(mfrow=c(3,2), mai = c(0.3, 0.3, 0.3, 0.3))
plot(pcadat$pc1,pcadat$pc2,main="PC1 PC2",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc1[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc2[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc2,pcadat$pc3,main="PC2 PC3",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc2[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc3[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc1,pcadat$pc3,main="PC1 PC3",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc1[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc3[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc1,pcadat$pc4,main="PC1 PC4",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc1[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc4[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc2,pcadat$pc4,main="PC2 PC4",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc2[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc4[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot.new()
legend("center", legend=DNA_plates, title="DNA extraction plate", fill=palette, cex=1.2, ncol=3)
dev.off()


# plot the phylogenetic diversity based on DNA plate (batch)

####### get sample size within each dna plate
cw_summary <- fpddat %>% 
  group_by(DNA_plate) %>% 
  tally()

#font size for pvalues 
your_font_size <- 3

# other fonts
My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 7),
  axis.title.y = element_text(size = 9))

batch_unroo <- ggboxplot(fpddat, x = "DNA_plate", y = "unrooted_pd",
                         color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 
batch_unroo
batch_bw <- ggboxplot(fpddat, x = "DNA_plate", y = "bwpd",
                      color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 
batch_bw

####### get sample size within each dna plate
cw_summary <- pcadat %>% 
  group_by(DNA_plate) %>% 
  tally()

batch_pc1 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc1",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-2.5, size = your_font_size) 
batch_pc1
batch_pc2 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc2",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=0, size = your_font_size) 
batch_pc2
batch_pc3 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc3",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.8, size = your_font_size) 
batch_pc3
batch_pc4 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc4",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.5, size = your_font_size) 
batch_pc4
batch_pc5 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=4, size = your_font_size)+
  ylim(4,7)
batch_pc5

# Extract the legend. Returns a gtable
for_legend_only <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                             color = "DNA_plate", palette = "jco")+
  guides(color=guide_legend(ncol=3)) +
  theme(legend.position=c(0.5,0.5),  
        plot.margin=unit(c(1,1,7,1),"lines")) +
  labs(fill="") 
for_legend_only
leg <- get_legend(for_legend_only)


figure <- grid.arrange(
  batch_unroo, batch_bw, batch_pc1, batch_pc2, batch_pc3, batch_pc4, batch_pc5, leg, nrow = 4, ncol = 2
)

pdf("out/batch_alpha_beta.pdf")
annotate_figure(figure,
                top = text_grob("Batch effect by alpha and beta diversity", color = "black", size = 14)
)
dev.off()


# adjusted p-values 

aov.out1 = aov(unrooted_pd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out1)
aov.out1 <- as.data.frame(res$DNA_plate)
aov.out1$type="unrooted_pd"

aov.out2 = aov(bwpd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out2)
aov.out2 <- as.data.frame(res$DNA_plate)
aov.out2$type="bwpd"

aov.out3 = aov(pc1 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out3)
aov.out3 <- as.data.frame(res$DNA_plate)
aov.out3$type="PC1"

aov.out4 = aov(pc2 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out4)
aov.out4 <- as.data.frame(res$DNA_plate)
aov.out4$type="PC2"

aov.out5 = aov(pc3 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out5)
aov.out5 <- as.data.frame(res$DNA_plate)
aov.out5$type="PC3"

aov.out6 = aov(pc4 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out6)
aov.out6 <- as.data.frame(res$DNA_plate)
aov.out6$type="PC4"

aov.out7 = aov(pc5 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out7)
aov.out7 <- as.data.frame(res$DNA_plate)
aov.out7$type="PC5"

aov.out1$type <- "unrooted_pd"
aov.out1$comparison <- rownames(aov.out1)
aov.out2$type <- "bwpd"
aov.out2$comparison <- rownames(aov.out2)
aov.out3$type <- "pc1"
aov.out3$comparison <- rownames(aov.out3)
aov.out4$type <- "pc2"
aov.out4$comparison <- rownames(aov.out4)
aov.out5$type <- "pc3"
aov.out5$comparison <- rownames(aov.out5)
aov.out6$type <- "pc4"
aov.out6$comparison <- rownames(aov.out6)
aov.out7$type <- "pc5"
aov.out7$comparison <- rownames(aov.out7)

all <- rbind(aov.out1,
      aov.out2,
      aov.out3,
      aov.out4,
      aov.out5,
      aov.out6,
      aov.out7)
all$padj_method <- "TukeyHSD"

addWorksheet(wb, "batch_pre_process")
writeData(wb, sheet = "batch_pre_process", all, rowNames = FALSE)

###########################################################################################

# Removes batch effect: 

PCA_well <- pcadat$DNA_well
PCA_batch <- pcadat$DNA_plate
head(pcadat)
pcadat<- data.matrix(pcadat[,3:7], rownames.force = NA)
pcadat_clean<-ComBat(dat=t(as.matrix(pcadat)),PCA_batch,mod=NULL)

PD_well <- fpddat$DNA_well
PD_batch <- fpddat$DNA_plate
fpddat<- data.matrix(fpddat[,3:7], rownames.force = NA)
fpddat_clean<-ComBat(dat=t(as.matrix(fpddat)),PD_batch,mod=NULL)

pcadat_clean <- t(pcadat_clean)
fpddat_clean <- t(fpddat_clean)

pcadat_clean <- as.data.frame(pcadat_clean)
fpddat_clean <- as.data.frame(fpddat_clean)

pcadat_clean <- unfactor(pcadat_clean[])
fpddat_clean <- unfactor(fpddat_clean[])

pcadat <- cbind(PCA_batch,PCA_well,pcadat_clean)
fpddat <- cbind(PD_batch,PD_well,fpddat_clean)

# rename cols to dna plate 
colnames(pcadat)[colnames(pcadat)=="PCA_batch"] <- "DNA_plate"
colnames(pcadat)[colnames(pcadat)=="PCA_well"] <- "DNA_well"
colnames(fpddat)[colnames(fpddat)=="PD_batch"] <- "DNA_plate"
colnames(fpddat)[colnames(fpddat)=="PD_well"] <- "DNA_well"


# save new un-batched data 
fwrite(x = pcadat, file = "out/pcadat_clean")
fwrite(x = fpddat, file = "out/fpddat_clean")


######### Plot batch effect after removal of batch effect: 

####### get sample size within each dna plate

#font size for pvalues 
your_font_size <- 3

# other fonts
My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 7),
  axis.title.y = element_text(size = 9))

####### get sample size within each dna plate
cw_summary <- fpddat %>% 
  group_by(DNA_plate) %>% 
  tally()

batch_unroo <- ggboxplot(fpddat, x = "DNA_plate", y = "unrooted_pd",
                         color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 
batch_unroo
batch_bw <- ggboxplot(fpddat, x = "DNA_plate", y = "bwpd",
                      color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 
batch_bw

cw_summary <- pcadat %>% 
  group_by(DNA_plate) %>% 
  tally()

batch_pc1 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc1",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-2.5, size = your_font_size) 
batch_pc1
batch_pc2 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc2",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=0, size = your_font_size) 
batch_pc2
batch_pc3 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc3",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.8, size = your_font_size) 
batch_pc3
batch_pc4 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc4",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.5, size = your_font_size) 
batch_pc4
batch_pc5 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=4, size = your_font_size)+
  ylim(4,7)
batch_pc5

# Extract the legend. Returns a gtable
for_legend_only <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                             color = "DNA_plate", palette = "jco")+
  guides(color=guide_legend(ncol=3)) +
  theme(legend.position=c(0.5,0.5),  
        plot.margin=unit(c(1,1,7,1),"lines")) +
  labs(fill="") 
for_legend_only
leg <- get_legend(for_legend_only)

figure <- grid.arrange(
  batch_unroo, batch_bw, batch_pc1, batch_pc2, batch_pc3, batch_pc4, batch_pc5, leg, nrow = 4, ncol = 2
)


# adjusted p-values 

aov.out1 = aov(unrooted_pd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out1)
aov.out1 <- as.data.frame(res$DNA_plate)
aov.out1$type="unrooted_pd"

aov.out2 = aov(bwpd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out2)
aov.out2 <- as.data.frame(res$DNA_plate)
aov.out2$type="bwpd"

aov.out3 = aov(pc1 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out3)
aov.out3 <- as.data.frame(res$DNA_plate)
aov.out3$type="PC1"

aov.out4 = aov(pc2 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out4)
aov.out4 <- as.data.frame(res$DNA_plate)
aov.out4$type="PC2"

aov.out5 = aov(pc3 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out5)
aov.out5 <- as.data.frame(res$DNA_plate)
aov.out5$type="PC3"

aov.out6 = aov(pc4 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out6)
aov.out6 <- as.data.frame(res$DNA_plate)
aov.out6$type="PC4"

aov.out7 = aov(pc5 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out7)
aov.out7 <- as.data.frame(res$DNA_plate)
aov.out7$type="PC5"

aov.out1$type <- "unrooted_pd"
aov.out1$comparison <- rownames(aov.out1)
aov.out2$type <- "bwpd"
aov.out2$comparison <- rownames(aov.out2)
aov.out3$type <- "pc1"
aov.out3$comparison <- rownames(aov.out3)
aov.out4$type <- "pc2"
aov.out4$comparison <- rownames(aov.out4)
aov.out5$type <- "pc3"
aov.out5$comparison <- rownames(aov.out5)
aov.out6$type <- "pc4"
aov.out6$comparison <- rownames(aov.out6)
aov.out7$type <- "pc5"
aov.out7$comparison <- rownames(aov.out7)

all <- rbind(aov.out1,
             aov.out2,
             aov.out3,
             aov.out4,
             aov.out5,
             aov.out6,
             aov.out7)
all$padj_method <- "TukeyHSD"

addWorksheet(wb, "batch_post_process")
writeData(wb, sheet = "batch_post_process", all, rowNames = FALSE)

pdf("out/NO_batch_alpha_beta.pdf")
annotate_figure(figure,
                top = text_grob("Batch effect after batch effect removal", color = "black", size = 14)
)
dev.off()


######################################################################################################

# 1.1 # merge alpha div data with metadata -> boggo
# 1.2 # merge beta div data with metadata -> coggo

# Merges with metadata and plots alpha diversity

# extracts the necessary columns from the metadata
mdat <- mdat %>%
  select(sample_name, isolation_source, collection_date, PigPen, Cohort, DNA_plate, DNA_well)


# alpha and beta div df formatting 
fpddat$DNA_plate <- gsub("P","plate_", fpddat$DNA_plate)
pcadat$DNA_plate <- gsub("P","plate_", pcadat$DNA_plate)

# Merge alpha and div to metadata one at a time

# merge metadata with alpha div 
boggo <- inner_join(fpddat,mdat)
NROW(boggo)
unique(boggo$Cohort)
sum(boggo$Cohort == "PosControl_ColiGuard")
sum(boggo$Cohort == "MockCommunity")
sum(boggo$Cohort == "PosControl_D-scour")
# NB:
# here cohorts "PosControl_ColiGuard" "PosControl_D-scour"   "MockCommunity" are present

# merge metadata with beta div 
coggo <- inner_join(pcadat,mdat)
NROW(coggo)
unique(coggo$Cohort)
# NB: 
# here cohorts "PosControl_ColiGuard" "PosControl_D-scour"   "MockCommunity" are missing
# and another 12 samples are missing (911-874=37 which is 8+9+8 pos controls and another 12 samples) 
# These are the non-matching samples: 
nomatch <- anti_join(boggo,coggo)
NROW(nomatch)
unique(nomatch$Cohort)
head(nomatch,37)

######################################################################################################

# 2   # plot samples distribution across plates (time and cohorts)


# use alpha div df (which includes the pos controls) to see how the samples are distributed in the plates 
# how are the samples distributed over the plates ? 
# some cohorts or dates over-represented in a plate? 

# frequency of DNA_plate by date

df1 <- setDT(boggo)[, .(Freq = .N), by = .(collection_date,DNA_plate)]

df1[order(df1$collection_date)]

p1 <- ggplot(df1, aes(fill=collection_date, y=Freq, x=DNA_plate)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "DNA extraction plate",
       y = "number of samples",
       fill = "collection date") +
  theme_bw()+
  theme(legend.position="top",
        axis.title.x=element_text(),
        legend.title=element_text(),
        axis.title.y=element_text())
p1

p2 <- ggplot(df1, aes(fill=DNA_plate, y=Freq, x=collection_date)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=palette)+
  labs(x = "sample collection date",
       y = "number of samples",
       fill = "DNA extraction plate") +
  theme_bw()+
  theme(legend.position="top",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())
p2

pdf("out/distribution_DNA_plate_time.pdf")
ggarrange(
  p1,p2,nrow=2, labels=c("A","B")
)
dev.off()

# frequency of date by Cohort

df1 <- setDT(boggo)[, .(Freq = .N), by = .(collection_date,Cohort)]

df1[order(df1$collection_date)]

pdf("out/distribution_cohorts_time.pdf")
ggplot(df1, aes(fill=Cohort, y=Freq, x=collection_date)) + 
  geom_bar(position="dodge", stat="identity")+
  labs(x = "time interval across trial",
       y = "number of samples",
       fill = "Cohort")+
  theme_bw()+
  theme(legend.position="top",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())
dev.off()


######################################################################################################

# 2.1 # boggo and coggo: formatting & averaging duplicates


# Now to get the stats and plot we need to remove (by averaging) 
# samples with identical collection date and identical isolation_source
# these are technical replicates

# aggregating dups for alpha
# select the necessary columns
boggo <- boggo %>%
  select(phylo_entropy,quadratic,unrooted_pd,rooted_pd,bwpd,isolation_source,
         collection_date, Cohort)
# necessary to remove and add later the pos and neg controls 
# otherwise when aggregating below we would end up with only 
# one replicate per pos/neg control
controls <- boggo %>%
  filter(Cohort == "MockCommunity"|
           Cohort == "PosControl_D-scour"|
           Cohort == "PosControl_ColiGuard"|
           Cohort == "NegativeControl")

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
boggo <- setDT(boggo)[, lapply(.SD, mean), by=c(names(boggo)[6:8]), .SDcols=cols]
# unite controls and samples 
boggo <- rbind(boggo, controls)

boggo$Cohort <- factor(boggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard",
                                "Mothers",
                                "MockCommunity",
                                "PosControl_D-scour",
                                "PosControlColiGuard"))

# aggregating dups for beta
# select the necessary columns
coggo <- coggo %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source,
         collection_date, Cohort)

# for beta we don't need to filter out then rbind the pos controls as we don't have them

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
coggo <- setDT(coggo)[, lapply(.SD, mean), by=c(names(coggo)[6:8]), .SDcols=cols]
coggo$Cohort <- factor(coggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard",
                                "Mothers",
                                "NegativeControl"))

######################################################################################################

# DELTAS 


# filtering out piglets that had dysentery
boggo1 <- boggo %>%
  filter(!isolation_source=="29665"|isolation_source=="29865"|isolation_source=="29702")

pigs_1 <- boggo1 %>%
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01") %>%
  select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_1)

hist(pigs_1$unrooted_pd,breaks=100)
hist(pigs_1$bwpd,breaks=100)

# unrooted_pd has low values indicative of spurious sample. 
# remove those rows
pigs_1 <- pigs_1 %>%
  filter(!unrooted_pd < 50) %>%
  filter(!bwpd >2.6) %>%
  filter(!bwpd <1.6)

###########################

pigs_2 <- boggo1 %>%
  filter(collection_date == "2017-02-06"|collection_date == "2017-02-07") %>%
  select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_2)

hist(pigs_2$unrooted_pd,breaks=100)
hist(pigs_2$bwpd,breaks=100)

# unrooted_pd has low values indicative of spurious sample. 
# remove those rows
pigs_2 <- pigs_2 %>%
  filter(!unrooted_pd < 90) %>%
  filter(!bwpd >2.3) %>%
  filter(!bwpd <1.8)

###########################

pigs_3 <- boggo1 %>%
  filter(collection_date == "2017-02-14") %>%
  select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_3)

hist(pigs_3$unrooted_pd)
hist(pigs_3$bwpd)

###########################

pigs_4 <- boggo1 %>%
  filter(collection_date == "2017-02-21") %>%
  select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_4)

hist(pigs_4$unrooted_pd,breaks=100)
hist(pigs_4$bwpd,100)

# unrooted_pd has low values indicative of spurious sample. 
# remove those rows
pigs_4 <- pigs_4 %>%
  filter(!unrooted_pd < 50) %>%
  filter(!bwpd < 1.6)

###########################

pigs_5 <- boggo1 %>%
  filter(collection_date == "2017-02-28") %>%
  select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_5)

hist(pigs_5$unrooted_pd)
hist(pigs_5$bwpd)


##############################################################################

# settings for plots: 

#font size for pvalues 
your_font_size <- 2 

My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

##############################################################################

# Ja31 vs Fe7

df1 <- merge(pigs_1,pigs_2, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

A_unroo <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_unroo, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))

A_bw <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_bw, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Ja31_vs_Fe7"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe7"
res2$type <- "bwpd"
A <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Ja31_vs_Fe7"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe7"
res2$type <- "bwpd"
A_adj <- rbind(res1,res2)


########## plots: 

cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

a1 <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylab("unrooted PD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

a2 <- ggboxplot(df1, x = "Cohort.x", y = "diff_bw", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylim(-35,15)+
  ylab("BWPD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

##############################################################################

# Fe7 vs Fe14

df1 <- merge(pigs_2,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

B_unroo <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_unroo, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))

B_bw <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_bw, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe7_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe14"
res2$type <- "bwpd"
B <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe7_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe14"
res2$type <- "bwpd"
B_adj <- rbind(res1,res2)


########## plots: 

cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

b1 <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylab("unrooted PD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

b2 <- ggboxplot(df1, x = "Cohort.x", y = "diff_bw", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylim(-35,25)+
  ylab("BWPD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

##############################################################################

# Fe14 vs Fe21

df1 <- merge(pigs_3,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

C_unroo <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_unroo, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))

C_bw <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_bw, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe14_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe21"
res2$type <- "bwpd"
C <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe14_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe21"
res2$type <- "bwpd"
C_adj <- rbind(res1,res2)


########## plots: 

cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

c1 <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylim(-75,30)+
  ylab("unrooted PD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

c2 <- ggboxplot(df1, x = "Cohort.x", y = "diff_bw", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylim(-23,25)+
  ylab("BWPD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

##############################################################################

# Fe21 vs Fe28

df1 <- merge(pigs_4,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

D_unroo <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_unroo, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))

D_bw <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_bw, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe21_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe21_vs_Fe28"
res2$type <- "bwpd"
D <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe21_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe21_vs_Fe28"
res2$type <- "bwpd"
D_adj <- rbind(res1,res2)


########## plots: 

cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

d1 <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylim(-100,50)+
  ylab("unrooted PD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

d2 <- ggboxplot(df1, x = "Cohort.x", y = "diff_bw", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylab("BWPD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

##############################################################################

# Ja31 vs Fe14

df1 <- merge(pigs_1,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

E_unroo <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_unroo, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))

E_bw <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_bw, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Ja31_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe14"
res2$type <- "bwpd"
E <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Ja31_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe14"
res2$type <- "bwpd"
E_adj <- rbind(res1,res2)


########## plots: 

cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

e1 <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylim(-75,30)+
  ylab("unrooted PD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

e2 <- ggboxplot(df1, x = "Cohort.x", y = "diff_bw", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylab("BWPD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

##############################################################################

# Fe7 vs Fe21

df1 <- merge(pigs_2,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

F_unroo <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_unroo, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))

F_bw <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_bw, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe7_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe21"
res2$type <- "bwpd"
FF <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe7_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe21"
res2$type <- "bwpd"
F_adj <- rbind(res1,res2)


########## plots: 

cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

f1 <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylab("unrooted PD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

f2 <- ggboxplot(df1, x = "Cohort.x", y = "diff_bw", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylab("BWPD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

##############################################################################

# Fe14 vs Fe28

df1 <- merge(pigs_3,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))

G_unroo <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_unroo, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))

G_bw <- df1 %>%
  group_by(Cohort.x) %>%
  summarize(Mean = mean(diff_bw, na.rm=TRUE),
            sd = sd(diff_unroo, na.rm=TRUE))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe14_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe28"
res2$type <- "bwpd"
G <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- as.data.frame(res1$p.value)
res2 <- as.data.frame(res2$p.value)
res1$time_delta <- "Fe14_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe28"
res2$type <- "bwpd"
G_adj <- rbind(res1,res2)


########## plots: 

cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

g1 <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylab("unrooted PD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


cw_summary <- df1 %>% 
  group_by(Cohort.x) %>% 
  tally()

g2 <- ggboxplot(df1, x = "Cohort.x", y = "diff_bw", color = "Cohort.x", 
                legend = "none") +
  My_Theme+
  ylim(-20,25)+
  ylab("BWPD - change (%)")+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 


##############################################################################

# this is for extracting the legend 
for_legend_only <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                             legend = "right")+
  scale_color_manual(labels = c("Control", 
                                "D-scour",
                                "ColiGuard",
                                "Neo",
                                "Neo+D",
                                "Neo+C"), 
                     values = c("#F8766D", 
                                "#B79F00",
                                "#00BA38",
                                "#00BFC4",
                                "#619CFF",
                                "#F564E3")) +
  guides(color=guide_legend("Cohort")) +
  My_Theme
leg <- get_legend(for_legend_only)

##############################################################################

# join the means and standard deviations ino a single dataframe : 

# unrooted_PD

colnames(A_unroo)[2:3] <- c("A-B_mean","A-B_sd")
colnames(E_unroo)[2:3] <- c("A-C_mean","A-C_sd")
colnames(B_unroo)[2:3] <- c("B-C_mean","B-C_sd")
colnames(F_unroo)[2:3] <- c("B-D_mean","B-D_sd")
colnames(C_unroo)[2:3] <- c("C-D_mean","C-D_sd")
colnames(G_unroo)[2:3] <- c("C-E_mean","C-E_sd")
colnames(D_unroo)[2:3] <- c("D-E_mean","D-E_sd")

unroo_means_and_sd <- join_all(list(A_unroo,
                                    E_unroo,
                                    B_unroo,
                                    F_unroo,
                                    C_unroo,
                                    G_unroo,
                                    D_unroo))

unroo_means_and_sd <- round_df(unroo_means_and_sd[,-1], 2)

unroo_means_and_sd$Cohort = c("Control",
                              "D-scour",
                              "ColiGuard",
                              "Neomycin",
                              "Neo+D-scour",
                              "Neo+ColiGuard")

unroo_means_and_sd$type="unrooted_PD"

# BWPD

colnames(A_bw)[2:3] <- c("A-B_mean","A-B_sd")
colnames(E_bw)[2:3] <- c("A-C_mean","A-C_sd")
colnames(B_bw)[2:3] <- c("B-C_mean","B-C_sd")
colnames(F_bw)[2:3] <- c("B-D_mean","B-D_sd")
colnames(C_bw)[2:3] <- c("C-D_mean","C-D_sd")
colnames(G_bw)[2:3] <- c("C-E_mean","C-E_sd")
colnames(D_bw)[2:3] <- c("D-E_mean","D-E_sd")

bwpd_means_and_sd <- join_all(list(A_bw,
                                   E_bw,
                                   B_bw,
                                   F_bw,
                                   C_bw,
                                   G_bw,
                                   D_bw))

bwpd_means_and_sd <- round_df(bwpd_means_and_sd[,-1], 2)

bwpd_means_and_sd$Cohort = c("Control",
                             "D-scour",
                             "ColiGuard",
                             "Neomycin",
                             "Neo+D-scour",
                             "Neo+ColiGuard")

bwpd_means_and_sd$type="BWPD"

# join unrooted_pd (means and sd) with bwpd (means and sd)
deltas_percent_change <- rbind(unroo_means_and_sd,bwpd_means_and_sd)

# add legend
deltas_percent_change$LEGEND <- c("NA","NA","NA",
                                  "A-B = Ja31 vs Fe7",
                                  "A-C = Ja31 vs Fe14",
                                  "B-C = Fe7 vs Fe14",
                                  "B-D = Fe7 vs Fe21",
                                  "C-D = Fe14 vs Fe21",
                                  "C-E = Fe14 vs Fe28",
                                  "D-E = Fe21 vs Fe28",
                                  "NA","NA")

# add data to workbook 
addWorksheet(wb, "deltas_percent_change")
writeData(wb, sheet = "deltas_percent_change", deltas_percent_change, rowNames = FALSE)


##############################################################################

# DELTAS p values  WITHOUT p-adjustment

# gather stats of the deltas: 
deltas_stats <- rbind(A,B,C,D,E,FF,G)
# convert rownames to first column
deltas_stats <- setDT(deltas_stats, keep.rownames = TRUE)[]
deltas_stats$test <- "pairwise t-test"
deltas_stats$padj_method <- "none"
# add data to workbook 
addWorksheet(wb, "deltas_p")
writeData(wb, sheet = "deltas_p", deltas_stats, rowNames = FALSE)

##############################################################################

# DELTAS p values WITH p-adjustment

# gather stats of the deltas: 
deltas_stats <- rbind(A_adj,B_adj,C_adj,D_adj,E_adj,F_adj,G_adj)
# convert rownames to first column
deltas_stats <- setDT(deltas_stats, keep.rownames = TRUE)[]
deltas_stats$test <- "pairwise t-test"
deltas_stats$padj_method <- "BH"
# add data to workbook 
addWorksheet(wb, "deltas_padj")
writeData(wb, sheet = "deltas_padj", deltas_stats, rowNames = FALSE)

##############################################################################

# all delta plots in two figures (unrooted and bwpd):

all_plots <- plot_grid(NULL, NULL, NULL, NULL, 
                       a1,e1,b1,f1,c1,g1,d1, leg, nrow = 3, 
                       labels = c("", "", "","", 
                                  "A-B", 
                                  "A-C", 
                                  "B-C", 
                                  "B-D", 
                                  "C-D", 
                                  "C-E", 
                                  "D-E", 
                                  ""),
                       ncol = 4)

pdf("out/cohorts_deltas_unrooted.pdf")
ggdraw() +
  draw_image(timeline_deltas_unroo, x = 0, y = 0.16) +
  draw_plot(all_plots)
dev.off()


all_plots <- plot_grid(NULL, NULL, NULL, NULL, 
                       a2,e2,b2,f2,c2,g2,d2, leg, nrow = 3, 
                       labels = c("", "", "","", 
                                  "A-B", 
                                  "A-C", 
                                  "B-C", 
                                  "B-D", 
                                  "C-D", 
                                  "C-E", 
                                  "D-E", 
                                  ""),
                       ncol = 4)

pdf("out/cohorts_deltas_bwpd.pdf")
ggdraw() +
  draw_image(timeline_deltas_bw, x = 0, y = 0.16) +
  draw_plot(all_plots)
dev.off()


######################################################################################################

# 3   # plot ALPHA diversity (all timepoints)

# ALPHA diversity overall (includes pos controls):

boggo<-inner_join(fpddat,mdat)
boggo <- boggo %>%
  filter(!isolation_source == "NegativeControl")

#pdf("out/alpha_phyloentropy.pdf",width=9,height=5)
#par(mar=(c(5, 10, 4, 2) +0.1))
#boxplot(boggo$phylo_entropy~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin+D-scour","Neomycin","Mothers","PosControl_ColiGuard","PosControl_D-scour","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Phylogenetic entropy",ylab=NULL,las=1)
#dev.off()

#pdf("out/alpha_unrooted.pdf",width=9,height=5)
#par(mar=(c(5, 10, 4, 2) +0.1))
#boxplot(boggo$unrooted_pd~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin+D-scour","Neomycin","Mothers","PosControl_ColiGuard","PosControl_D-scour","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Unrooted PD",ylab=NULL,las=1)
#dev.off()

#pdf("out/alpha_bwpd.pdf",width=9,height=5)
#par(mar=(c(5, 10, 4, 2) +0.1))
#boxplot(boggo$bwpd~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin+D-scour","Neomycin","Mothers","PosControl_ColiGuard","PosControl_D-scour","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Balance-weighted PD",ylab=NULL,las=1)
#dev.off()


# boxplots again for alpha diversity, to be plotted in the same pdf
p1 <- ggplot(boggo, aes(x=fct_inorder(Cohort), y=phylo_entropy)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()
p2 <- ggplot(boggo, aes(x=fct_inorder(Cohort), y=bwpd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()
p3 <- ggplot(boggo, aes(x=fct_inorder(Cohort), y=unrooted_pd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()
p3

pdf("out/alpha_BWPD&unrooted.pdf")
plot_grid(p2,p3, align = "hv", nrow=2, labels = "auto")
dev.off()

pdf("out/alpha_all.pdf")
plot_grid(p1,p2,p3, align = "hv", nrow=3, labels = "auto")
dev.off()

##############################################

a <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Mothers") %>% 
  dplyr::summarise(min = min(unrooted_pd)
            ,max = max(unrooted_pd)
            ,mean = mean(unrooted_pd)
            ,n = n()
            ,sd = sd(unrooted_pd)
            ,q25 = quantile(unrooted_pd, .25)
            ,q75 = quantile(unrooted_pd, .75)) 

b <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Control") %>% 
  filter(!Cohort=="D-scour") %>% 
  filter(!Cohort=="ColiGuard") %>% 
  filter(!Cohort=="Neomycin") %>% 
  filter(!Cohort=="Neomycin+D-scour") %>% 
  filter(!Cohort=="Neomycin+ColiGuard") %>% 
  dplyr::summarise(min = min(unrooted_pd)
            ,max = max(unrooted_pd)
            ,mean = mean(unrooted_pd)
            ,n = n()
            ,sd = sd(unrooted_pd)
            ,q25 = quantile(unrooted_pd, .25)
            ,q75 = quantile(unrooted_pd, .75)) 

c <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Mothers") %>% 
  dplyr::summarise(min = min(bwpd)
            ,max = max(bwpd)
            ,mean = mean(bwpd)
            ,n = n()
            ,sd = sd(bwpd)
            ,q25 = quantile(bwpd, .25)
            ,q75 = quantile(bwpd, .75)) 

d <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Control") %>% 
  filter(!Cohort=="D-scour") %>% 
  filter(!Cohort=="ColiGuard") %>% 
  filter(!Cohort=="Neomycin") %>% 
  filter(!Cohort=="Neomycin+D-scour") %>% 
  filter(!Cohort=="Neomycin+ColiGuard") %>% 
  dplyr::summarise(min = min(bwpd)
            ,max = max(bwpd)
            ,mean = mean(bwpd)
            ,n = n()
            ,sd = sd(bwpd)
            ,q25 = quantile(bwpd, .25)
            ,q75 = quantile(bwpd, .75)) 

e <- boggo %>% group_by(Cohort) %>% 
  dplyr::summarise(min = min(unrooted_pd)
            ,max = max(unrooted_pd)
            ,mean = mean(unrooted_pd)
            ,n = n()
            ,sd = sd(unrooted_pd)
            ,q25 = quantile(unrooted_pd, .25)
            ,q75 = quantile(unrooted_pd, .75))

f <- boggo %>% group_by(Cohort) %>% 
  dplyr::summarise(min = min(bwpd)
            ,max = max(bwpd)
            ,mean = mean(bwpd)
            ,n = n()
            ,sd = sd(bwpd)
            ,q25 = quantile(bwpd, .25)
            ,q75 = quantile(bwpd, .75)) 

a$Cohort="piglets"
b$Cohort="mothers"
c$Cohort="piglets"
d$Cohort="mothers"
e$Cohort="all"
f$Cohort="all"

a$type="unrooted"
b$type="unrooted"
c$type="bwpd"
d$type="bwpd"
e$type="unrooted"
f$type="bwpd"

means <- rbind(a,b,c,d,e,f)
means$collection_date = "all"
# these are added later below into a df 

##############################################

# alpha_timeseries_all_cohorts

doggo <- boggo
doggo$collection_date <- as.character(doggo$collection_date)

doggo <- doggo %>% filter(collection_date == "2017-01-31" |
    collection_date == "2017-02-07" |
    collection_date == "2017-02-14" |
    collection_date == "2017-02-21" )

# reordering
doggo$Cohort <- factor(doggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"))
# reordering
doggo$collection_date <- factor(doggo$collection_date,
                                levels=c("2017-01-31",
                                         "2017-02-07",
                                         "2017-02-14",
                                         "2017-02-21" ))

my_comparisons = list( c("2017-01-31", "2017-02-07"), 
                        c("2017-02-07", "2017-02-14"), 
                        c("2017-02-14", "2017-02-21"),
                        c("2017-02-07", "2017-02-21"))

# general time change - unrooted
summs_unroo <- doggo %>% group_by(collection_date,Cohort) %>% 
  dplyr::summarise(min = min(unrooted_pd)
                   ,max = max(unrooted_pd)
            ,mean = mean(unrooted_pd)
            ,n = n()
            ,sd = sd(unrooted_pd)
            ,q25 = quantile(unrooted_pd, .25)
            ,q75 = quantile(unrooted_pd, .75)) 

gen_unrooted <- ggplot(summs_unroo, aes(x=collection_date, y=mean, group=Cohort, color=Cohort)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1,
                position=position_dodge(0.5)) +
  geom_line() + geom_point()+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = "collection date",
       y = "unrooted PD - mean")
gen_unrooted


# general time change - unrooted
summs_bw <- doggo %>% group_by(collection_date,Cohort) %>% 
  dplyr::summarise(min = min(bwpd)
            ,max = max(bwpd)
            ,mean = mean(bwpd)
            ,n = n()
            ,sd = sd(bwpd)
            ,q25 = quantile(bwpd, .25)
            ,q75 = quantile(bwpd, .75)) 

gen_bwpd <- ggplot(summs_bw, aes(x=collection_date, y=mean, group=Cohort, color=Cohort)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1,
                position=position_dodge(0.5)) +
  geom_line() + geom_point()+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = "collection date",
       y = "BWPD - mean")
gen_bwpd


pdf("out/time_alpha.pdf")
grid.arrange(
  gen_unrooted, gen_bwpd, nrow = 2
)
dev.off()

# unrooted pd - fill: collection date
#pdf("out/cohorts_unrooted.pdf",width=9,height=5)
# unrooted_time <- ggplot(doggo, aes(x=Cohort, y=unrooted_pd, 
#                                    fill=collection_date)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_boxplot() +
#   scale_fill_discrete(name = "collection date") + 
#   scale_y_continuous(limits = c(60, 170))
# unrooted_time
#dev.off()

# bwpd - fill: collection date
#pdf("out/alpha_bwpd_fill_cohorts.pdf",width=9,height=5)
# bwpd_time <- ggplot(doggo, aes(x=Cohort, y=bwpd, 
#                                fill=collection_date)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_boxplot() +
#   scale_fill_discrete(name = "collection date") + 
#   scale_y_continuous(limits = c(1.75, 2.8))
# bwpd_time
#dev.off()

#pdf("out/cohorts_unrooted&bwpd.pdf")
#grid.arrange(
#  unrooted_time, bwpd_time, nrow = 2,
#  top = "Alpha diversity"
#)
#dev.off()

#########################################################################

# stats 

unroo <- doggo %>%
  group_by(collection_date,Cohort) %>%
  dplyr::summarise(min = min(unrooted_pd)
            ,max = max(unrooted_pd)
            ,mean = mean(unrooted_pd)
            ,sd = sd(unrooted_pd)
            ,n = n()
            ,q25 = quantile(unrooted_pd, .25)
            ,q75 = quantile(unrooted_pd, .75)) 

bw <- doggo %>%
  group_by(collection_date,Cohort) %>%
  dplyr::summarise(min = min(bwpd)
            ,max = max(bwpd)
            ,mean = mean(bwpd)
            ,sd = sd(bwpd)
            ,n = n()
            ,q25 = quantile(bwpd, .25)
            ,q75 = quantile(bwpd, .75)) 

# add data to workbook 

bw$type="bwpd"
unroo$type="unrooted_pd"
bw <- as.data.frame(bw)
unroo <- as.data.frame(unroo)
both <- rbind(bw,unroo,means)
addWorksheet(wb, "alpha_means")
writeData(wb, sheet = "alpha_means", both, rowNames = FALSE)

################


stat.test_unroo <- doggo %>%
  group_by(Cohort) %>%
  t_test(unrooted_pd ~ collection_date) %>%
  adjust_pvalue(method="BH") %>%
  mutate(y.position=rep(seq(150,250,length.out=6),6)) %>%
  mutate_if(is.numeric, round, digits = 4)


unroo <- ggboxplot(doggo, x = "collection_date", y = "unrooted_pd",
               color = "collection_date", palette = "jco",
               add = "jitter",facet.by = "Cohort", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(0,250)+
  stat_pvalue_manual(stat.test_unroo, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = 3)

#########

stat.test_bwpd <- doggo %>%
  group_by(Cohort) %>%
  t_test(bwpd ~ collection_date) %>%
  adjust_pvalue(method="BH") %>%
  mutate(y.position=rep(seq(2.5,2.8,length.out=6),6)) %>%
  mutate_if(is.numeric, round, digits = 4)

bw <- ggboxplot(doggo, x = "collection_date", y = "bwpd",
          color = "collection_date", palette = "jco",
          add = "jitter",facet.by = "Cohort", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(1.3,3)+
  stat_pvalue_manual(stat.test_bwpd, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = 3)


pdf("out/time_alpha_cohorts.pdf")
ggarrange(unroo,bw,nrow=2,labels=c("A","B"))
dev.off()


# add data to workbook 

both <- rbind(stat.test_bwpd,
      stat.test_unroo)
both$padj_method <- "BH"

addWorksheet(wb, "alpha_cohorts")
writeData(wb, sheet = "alpha_cohorts", both, rowNames = FALSE)

######################################################################################################

# 4   # plot BETA diversity (all timepoints)

# BETA diversity overall (does not include pos controls):

# plot

pdf("out/pca.pdf")
plot(coggo$pc1,coggo$pc3,main="phylosift edge PCA on pigs",xlab="PC1",ylab="PC3",pch=NA,type="p")
cohorts=unique(sort(coggo$Cohort))
for(coho in 1:length(cohorts)){
  points(coggo$pc1[coggo$Cohort==cohorts[coho]],coggo$pc3[coggo$Cohort==cohorts[coho]],col=coho,pch=1)
}
legend("topleft",legend=cohorts, fill=palette(), cex=0.8, ncol=2)
dev.off()

#2*ylen/3
color_legend <- function(x, y, xlen, ylen, main, tiks, colors){
  text(x, y+.6, main, adj=c(0,0), cex=1.3)
  color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.8)
}

rbow <- rainbow(40, end=0.7, alpha=0.7)

pdf("out/time_beta.pdf")
plot(coggo$pc1[coggo$Cohort!="Mothers"&coggo$Cohort!="NegativeControl"],coggo$pc2[coggo$Cohort!="Mothers"&coggo$Cohort!="NegativeControl"],main="beta diversity (phylosift edge PCA)",xlab="PC1",ylab="PC2",type="p",col=rbow[as.Date(coggo$collection_date[coggo$Cohort!="Mothers"&coggo$Cohort!="NegativeControl"])-as.Date("2017-01-29 00:00:00")])
legvec <- c(0,10,20,30,40)
color_legend( -2.9, 4, 3.5, 1.5, "trial days:", legvec, rbow)
dev.off()

#################################

# timeseries within cohort 

pdf("out/time_beta_cohorts.pdf") 
par(mfrow=c(3,2), mai = c(0.4, 0.4, 0.4, 0.4))
plot(coggo$pc1[coggo$Cohort=="Control"],
     coggo$pc2[coggo$Cohort=="Control"],
     main="PC1 PC2 Control",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc3[coggo$Cohort=="Control"],
     coggo$pc4[coggo$Cohort=="Control"],
     main="PC3 PC4 Control",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc1[coggo$Cohort=="D-scour"],
     coggo$pc2[coggo$Cohort=="D-scour"],
     main="PC1 PC2 D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc3[coggo$Cohort=="D-scour"],
     coggo$pc4[coggo$Cohort=="D-scour"],
     main="PC3 PC4 D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc1[coggo$Cohort=="ColiGuard"],
     coggo$pc2[coggo$Cohort=="ColiGuard"],
     main="PC1 PC2 ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc3[coggo$Cohort=="ColiGuard"],
     coggo$pc4[coggo$Cohort=="ColiGuard"],
     main="PC3 PC4 ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
par(mfrow=c(3,2), mai = c(0.4, 0.4, 0.4, 0.4))
plot(coggo$pc1[coggo$Cohort=="Neomycin"],
     coggo$pc2[coggo$Cohort=="Neomycin"],
     main="PC1 PC2 Neomycin",
     xlab="PC1",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc3[coggo$Cohort=="Neomycin"],
     coggo$pc4[coggo$Cohort=="Neomycin"],
     main="PC3 PC4 Neomycin",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc1[coggo$Cohort=="Neomycin+D-scour"],
     coggo$pc2[coggo$Cohort=="Neomycin+D-scour"],
     main="PC1 PC2 Neomycin+D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Neomycin+D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc3[coggo$Cohort=="Neomycin+D-scour"],
     coggo$pc4[coggo$Cohort=="Neomycin+D-scour"],
     main="PC3 PC4 Neomycin+D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Neomycin+D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc1[coggo$Cohort=="Neomycin+ColiGuard"],
     coggo$pc2[coggo$Cohort=="Neomycin+ColiGuard"],
     main="PC1 PC2 Neomycin+ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
plot(coggo$pc3[coggo$Cohort=="Neomycin+ColiGuard"],
     coggo$pc4[coggo$Cohort=="Neomycin+ColiGuard"],
     main="PC3 PC4 Neomycin+ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(coggo$collection_date
                               [coggo$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
dev.off()


######################################################################################################

# 5   # plot ALPHA diversity (at pig trial start) 

# ALPHA diversity in piglets at the start of the trial 

startDF <- boggo 

startDF <- startDF %>% filter(
  collection_date == "2017-01-31" |
    collection_date == "2017-02-01")  %>%
  select(phylo_entropy,quadratic,unrooted_pd,rooted_pd,bwpd,isolation_source)

# as we have 160 samples for 126 piglets at the startof the trial. this is because we have duplicates
# samples have been taken from the same animals twice
length(startDF$isolation_source)
length(unique(startDF$isolation_source))
# we need to average 
head(startDF)
# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
startDF <- setDT(startDF)[, lapply(.SD, mean), by=c(names(startDF)[6]), .SDcols=cols]
# now we have the right number: 
length(startDF$isolation_source)
length(unique(startDF$isolation_source))


#########################

startDF1 <- merge(startDF,details, by.x="isolation_source",by.y="pig")

startDF1 <- startDF1 %>%
  select(phylo_entropy,unrooted_pd,bwpd,isolation_source,nurse,stig)


# plots

cw_summary <- startDF1 %>% 
  group_by(nurse) %>% 
  tally()

a <- ggplot(startDF1, aes(x=nurse, y=bwpd, group=nurse)) + 
  labs(title = "Piglets alpha diversity (BWPD)",
       subtitle = "Grouped by nurse mother") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(nurse, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11))

b <- ggplot(startDF1, aes(x=nurse, y=unrooted_pd, group=nurse)) + 
  labs(title = "Piglets alpha diversity (unrooted)",
       subtitle = "Grouped by nurse mother") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(nurse, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11)) +
  scale_y_continuous(limits=c(60,160))

#pdf("out/nurse_alpha.pdf")
# ggarrange(a, b, 
#           labels = c("A", "B"),
#           ncol = 1, nrow = 2)
#dev.off()


############

cw_summary <- startDF1 %>% 
  group_by(stig) %>% 
  tally()

c <- ggplot(startDF1, aes(x=stig, y=bwpd, group=stig)) + 
  labs(title = "Piglets alpha diversity (BWPD)",
       subtitle = "Grouped by stig mother") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(stig, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11))

d <- ggplot(startDF1, aes(x=stig, y=unrooted_pd, group=stig)) + 
  labs(title = "Piglets alpha diversity (unrooted)",
       subtitle = "Grouped by stig mother") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(stig, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11)) +
  scale_y_continuous(limits=c(60,160))


#pdf("out/stig_alpha.pdf")
# ggarrange(c, d, 
#           labels = c("A", "B"),
#           ncol = 1, nrow = 2)
#dev.off()


##################

# same but putting nurses and stigs on the same plot, dividing BWPD from unrooted
pdf("out/nurse&stig_BWPD.pdf")
ggarrange(c, a, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

pdf("out/nurse&stig_unrooted.pdf")
ggarrange(d, b, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

##############################################

# do piglets at arrival cluster by unrooted and bwpd
# based on breed/line/birth day? 

startDF2 <- merge(startDF,details, by.x="isolation_source",by.y="pig")

startDF2$BIRTH_DAY <- as.character(startDF2$BIRTH_DAY)
startDF2$BIRTH_DAY <- factor(startDF2$BIRTH_DAY, 
                             levels=c("2017-01-06", 
                                      "2017-01-07", 
                                      "2017-01-08",
                                      "2017-01-09",
                                      "2017-01-10",
                                      "2017-01-11"))
startDF2$LINE <- as.character(startDF2$LINE)

# plots

# by breed

####### get sample size within each breed group:

cw_summary <- startDF2 %>% 
  group_by(breed) %>% 
  tally()

# breed - unrooted 
breed_unrooted_plot <- ggboxplot(startDF2, x = "breed", y = "unrooted_pd",
                                 color = "breed", palette = "jco",
                                 add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(50,160)+
  geom_text(data = cw_summary,
            aes(breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# breed - bwpd

breed_bwpd_plot <- ggboxplot(startDF2, x = "breed", y = "bwpd",
                             color = "breed", palette = "jco",
                             add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(),legend.position="none")+
  ylim(1.5,2.7)+
  geom_text(data = cw_summary,
            aes(breed, Inf, label = n), vjust="inward")+
  stat_compare_means(method = "kruskal.test", label.y=1.5) 

pdf("out/breed_alpha.pdf")
grid.arrange(
  breed_unrooted_plot, breed_bwpd_plot, nrow = 2
)
dev.off()

# by birth day

####### get sample size within each bday group:

cw_summary <- startDF2 %>% 
  group_by(BIRTH_DAY) %>% 
  tally()

# bday - unrooted 
bday_unrooted_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "unrooted_pd",
                                color = "BIRTH_DAY", palette = "jco",
                                add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(50,160)+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# bday - bwpd 
bday_bwpd_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "bwpd",
                            color = "BIRTH_DAY", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  ylim(1.5,2.7)+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=1.5)  # Add pairwise comparisons p-value

pdf("out/bday_alpha.pdf")
grid.arrange(
  bday_unrooted_plot, bday_bwpd_plot, nrow = 2
)
dev.off()


# by line

####### get sample size within each bday group:

cw_summary <- startDF2 %>% 
  group_by(LINE) %>% 
  tally()

# line - unrooted 
LINE_unrooted_plot <- ggboxplot(startDF2, x = "LINE", y = "unrooted_pd",
                                color = "LINE", palette = "jco",
                                add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(50,160)+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# line - bwpd 
LINE_bwpd_plot <- ggboxplot(startDF2, x = "LINE", y = "bwpd",
                            color = "LINE", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  ylim(1.5,2.75)+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=1.5)  # Add pairwise comparisons p-value

pdf("out/line_alpha.pdf")
grid.arrange(
  LINE_unrooted_plot, LINE_bwpd_plot, nrow = 2
)
dev.off()

##################

# bday AND breed:
# with age unrooted pd decreases and bwpd increases 
# conclusion comes from two breeds as each bday is not represented 
# equally by the 4 breeds 
my_comparisons <- list(  c("2017-01-07", "2017-01-09"),  
                         c("2017-01-08", "2017-01-10"), 
                         c("2017-01-09", "2017-01-11"), 
                         c("2017-01-10", "2017-01-11") )

startDF2_sub <- startDF2 %>%
  filter(!breed=="Landrace x Cross bred (LW x D)")
startDF2_sub <- startDF2_sub %>%
  filter(!breed=="Large white x Duroc")


p1 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "unrooted_pd",
               color = "BIRTH_DAY", palette = "jco",
               add = "jitter",
               facet.by = "breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(50,230)+
  stat_compare_means(comparisons = my_comparisons)

p2 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "bwpd",
               color = "BIRTH_DAY", palette = "jco",
               add = "jitter",
               facet.by = "breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(1.6,3.4)+
  stat_compare_means(comparisons = my_comparisons)

pdf("out/bday_bybreed_alpha.pdf",width=9,height=5)
ggarrange(p1,p2,nrow=2)
dev.off()


######################################################################################################

# 6   # plot BETA diversity (at pig trial start) 

# BETA diversity in piglets at the start of the trial 

startDF <- coggo %>% filter(
  collection_date == "2017-01-31" |
    collection_date == "2017-02-01")  %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source)

# as we have 160 samples for 126 piglets at the startof the trial. this is because we have duplicates
# samples have been taken from the same animals twice
length(startDF$isolation_source)
length(unique(startDF$isolation_source))
# we need to average 
head(startDF)
# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
startDF <- setDT(startDF)[, lapply(.SD, mean), by=c(names(startDF)[6]), .SDcols=cols]
# now we have the right number: 
length(startDF$isolation_source)
length(unique(startDF$isolation_source))


#########################

startDF1 <- merge(startDF,details, by.x="isolation_source",by.y="pig")

startDF1 <- startDF1 %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source,nurse,stig)

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.title.x=element_text(colour="black",size=8),
             axis.title.y=element_text(colour="black",size=8),
             axis.text.x=element_text(colour="black",size=8),
             axis.text.y=element_text(colour="black",size=8),
             axis.ticks=element_line(colour="black"),
             legend.position="right",
             legend.text=element_text(size=8),
             legend.title=element_text(size=8),
             plot.margin=unit(c(0.3,0.3,0.3,0.3),"line"))

# nurses

startDF1_unique <- startDF1 %>% group_by(isolation_source) %>% slice(1)

# remove rows were nurse is unique (can't be plotted a pca)
startDF1_unique <- subset(startDF1_unique,duplicated(nurse) | duplicated(nurse, fromLast=TRUE))

# order alphabetically by nurse
startDF1_unique <- startDF1_unique[order(startDF1_unique$nurse),]

rownames(startDF1_unique) <- 1:nrow(startDF1_unique)
length(unique(startDF1_unique$nurse))

# subsetting to have 5 nurses in each 
startDF11 <- startDF1_unique[1:27,]
startDF12 <- startDF1_unique[28:47,]
startDF13 <- startDF1_unique[48:63,]
startDF14 <- startDF1_unique[64:80,]
startDF15 <- startDF1_unique[81:98,]
startDF16 <- startDF1_unique[99:122,]

p11<-ggplot(startDF11,aes(x=pc1,y=pc2,color=nurse ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p12<-ggplot(startDF12,aes(x=pc1,y=pc2,color=nurse ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p13<-ggplot(startDF13,aes(x=pc1,y=pc2,color=nurse ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p14<-ggplot(startDF14,aes(x=pc1,y=pc2,color=nurse ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p15<-ggplot(startDF15,aes(x=pc1,y=pc2,color=nurse ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p16<-ggplot(startDF16,aes(x=pc1,y=pc2,color=nurse ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, 
               level = 0.80)
p16

pdf("out/nurse_PC1PC2.pdf")
grid.arrange(p11,p12,p13,p14,p15,p16, nrow=3,ncol=2)
dev.off()

######################

# stigs 

startDF1_unique <- startDF1 %>% group_by(isolation_source) %>% slice(1)

# remove rows were nurse is unique (can't be plotted a pca)
startDF1_unique <- subset(startDF1_unique,duplicated(stig) | duplicated(stig, fromLast=TRUE))

# order alphabetically by nurse
startDF1_unique <- startDF1_unique[order(startDF1_unique$stig),]
rownames(startDF1_unique) <- 1:nrow(startDF1_unique)
length(unique(startDF1_unique$stig))

# subsetting to have 5 stigs in each 
startDF11 <- startDF1_unique[1:30,]
startDF12 <- startDF1_unique[31:58,]
startDF13 <- startDF1_unique[59:72,]
startDF14 <- startDF1_unique[73:85,]
startDF15 <- startDF1_unique[86:108,]
startDF16 <- startDF1_unique[109:123,]

p11<-ggplot(startDF11,aes(x=pc1,y=pc2,color=stig ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p12<-ggplot(startDF12,aes(x=pc1,y=pc2,color=stig ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p13<-ggplot(startDF13,aes(x=pc1,y=pc2,color=stig ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p14<-ggplot(startDF14,aes(x=pc1,y=pc2,color=stig ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p15<-ggplot(startDF15,aes(x=pc1,y=pc2,color=stig ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p16<-ggplot(startDF16,aes(x=pc1,y=pc2,color=stig ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, 
               level = 0.80)
p16

pdf("out/stig_PC1PC2.pdf")
grid.arrange(p11,p12,p13,p14,p15,p16, nrow=3,ncol=2)
dev.off()

##############################################

# do piglets at arrival cluster by PCA
# based on breed/line/birth day? 

startDF2 <- merge(startDF,details, by.x="isolation_source",by.y="pig")

startDF2$BIRTH_DAY <- as.character(startDF2$BIRTH_DAY)
startDF2$BIRTH_DAY <- factor(startDF2$BIRTH_DAY, 
                             levels=c("2017-01-06", 
                                      "2017-01-07", 
                                      "2017-01-08",
                                      "2017-01-09",
                                      "2017-01-10",
                                      "2017-01-11"))

startDF2$LINE <- as.character(startDF2$LINE)


# by breed

####### get sample size within each breed group:

cw_summary <- startDF2 %>% 
  group_by(breed) %>% 
  tally()

# breed - PCA
breed_PC1_plot <- ggboxplot(startDF2, x = "breed", y = "pc1",
                            color = "breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.5)  # Add pairwise comparisons p-value
breed_PC1_plot
breed_PC2_plot <- ggboxplot(startDF2, x = "breed", y = "pc2",
                            color = "breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=1)  # Add pairwise comparisons p-value
breed_PC2_plot
breed_PC3_plot <- ggboxplot(startDF2, x = "breed", y = "pc3",
                            color = "breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.5)  # Add pairwise comparisons p-value
breed_PC3_plot
breed_PC4_plot <- ggboxplot(startDF2, x = "breed", y = "pc4",
                            color = "breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.2)  # Add pairwise comparisons p-value
breed_PC4_plot
breed_PC5_plot <- ggboxplot(startDF2, x = "breed", y = "pc5",
                            color = "breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(breed, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value
breed_PC5_plot

pdf("out/breed_beta.pdf")
grid.arrange(
  breed_PC1_plot, breed_PC2_plot, breed_PC3_plot, breed_PC4_plot, breed_PC5_plot, nrow = 3, ncol=2
)
dev.off()

# by birth day

####### get sample size within each bday group:

cw_summary <- startDF2 %>% 
  group_by(BIRTH_DAY) %>% 
  tally()

# bday - PCA
bday_PC1_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc1",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.2)  # Add pairwise comparisons p-value
bday_PC1_plot
bday_PC2_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc2",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=1)  # Add pairwise comparisons p-value
bday_PC2_plot
bday_PC3_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc3",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.5)  # Add pairwise comparisons p-value
bday_PC3_plot
bday_PC4_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc4",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1)  # Add pairwise comparisons p-value
bday_PC4_plot
bday_PC5_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc5",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value
bday_PC5_plot

pdf("out/bday_beta.pdf")
grid.arrange(
  bday_PC1_plot, bday_PC2_plot, bday_PC3_plot, bday_PC4_plot, bday_PC5_plot, nrow = 3, ncol=2
)
dev.off()


# by line

####### get sample size within each bday group:

cw_summary <- startDF2 %>% 
  group_by(LINE) %>% 
  tally()

# line - PCA
line_PC1_plot <- ggboxplot(startDF2, x = "LINE", y = "pc1",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.2)  # Add pairwise comparisons p-value
line_PC1_plot
line_PC2_plot <- ggboxplot(startDF2, x = "LINE", y = "pc2",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=.2)  # Add pairwise comparisons p-value
line_PC2_plot
line_PC3_plot <- ggboxplot(startDF2, x = "LINE", y = "pc3",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.4)  # Add pairwise comparisons p-value
line_PC3_plot
line_PC4_plot <- ggboxplot(startDF2, x = "LINE", y = "pc4",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1)  # Add pairwise comparisons p-value
line_PC4_plot
line_PC5_plot <- ggboxplot(startDF2, x = "LINE", y = "pc5",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value
line_PC5_plot

pdf("out/line_beta.pdf")
grid.arrange(
  line_PC1_plot, line_PC2_plot, line_PC3_plot, line_PC4_plot, line_PC5_plot, nrow = 3, ncol=2
)
dev.off()

##################

# bday AND breed:

my_comparisons <- list(  c("2017-01-07", "2017-01-09"),  
                         c("2017-01-08", "2017-01-10"), 
                         c("2017-01-09", "2017-01-11"), 
                         c("2017-01-10", "2017-01-11") )

# need to exclude two breeds as these two breeds don't have enough
# age groups to be plotted and compared 
startDF2_sub <- startDF2 %>%
  filter(!breed=="Landrace x Cross bred (LW x D)")
startDF2_sub <- startDF2_sub %>%
  filter(!breed=="Large white x Duroc")


p1 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc1",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(0,4)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p2 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc2",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(0,7)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p3 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc3",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(-2,3)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p4 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc4",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(-1,2.2)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p5 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc5",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(4.7,6.6)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)

pdf("out/bday_bybreed_beta.pdf")
grid.arrange(
  p1,p2, nrow = 2, ncol = 1
)
grid.arrange(
  p3,p4, nrow = 2, ncol = 1
)
grid.arrange(
  p5, nrow = 2, ncol = 1
)
dev.off()


######################################################################################################

# distribution of breeds, birth days across cohorts

finalDF <- inner_join(boggo,coggo)
NROW(finalDF)

df <- merge(finalDF,details, by.x="isolation_source",by.y="pig")
head(df)

# distribution of breeds and bdays across cohorts
df1 <- setDT(df)[, .(Freq = .N), by = .(BIRTH_DAY,breed,Cohort)]
df1[order(df1$breed)]

p1 <- ggplot(df1, aes(fill=BIRTH_DAY, y=Freq, x=Cohort)) + 
  geom_bar(position="stack", stat="identity")+
  theme(legend.position="none")+
  facet_wrap(.~breed,scales="free")
p1

# distribution of Cohorts across bdays
df1 <- setDT(df)[, .(Freq = .N), by = .(Cohort,BIRTH_DAY)]
df1[order(df1$Cohort)]

p2 <- ggplot(df1, aes(fill=sort(BIRTH_DAY), y=Freq, x=Cohort)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "Cohort",
       y = "number of samples",
       fill = "birth day") +
  theme_bw()+
  theme(legend.position="right",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())
p2

pdf("out/distribution_breeds&bday_cohorts.pdf")
ggarrange(p1,p2,nrow=2,labels=c("A","B"))
dev.off()

# distribution of Cohorts across DNA plates
df1 <- setDT(df)[, .(Freq = .N), by = .(Cohort,DNA_plate)]
df1[order(df1$Cohort)]

pdf("out/distribution_cohorts_DNA_plate.pdf")
p2 <- ggplot(df1, aes(fill=DNA_plate, y=Freq, x=Cohort)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=palette)+
  labs(x = "Cohort",
       y = "number of samples",
       fill = "DNA extraction plate") +
  theme_bw()+
  theme(legend.position="top",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())
p2
dev.off()


######################################################################################################

# 8   # p-values

# PVALUES - 


# merge alpha&beta div (finalDF) to details and details metadata

df <- merge(finalDF,details, by.x="isolation_source",by.y="pig")
head(df)

# rename collection dates to time intervals 
# iM <- "2017-01-30"
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-01-30", 
  replacement = "iM", 
  fixed = TRUE)
# i0 <- "2017-01-31" "2017-02-01" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-01-31", 
  replacement = "i1", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-01", 
  replacement = "i1", 
  fixed = TRUE)

# i1 <- "2017-02-03" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-03", 
  replacement = "i2.1", 
  fixed = TRUE)

# i2 <- "2017-02-06" "2017-02-07" "2017-02-08"
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-06", 
  replacement = "i2.2", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-07", 
  replacement = "i2.2", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-08", 
  replacement = "i2.2", 
  fixed = TRUE)

# i3 <- "2017-02-10" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-10", 
  replacement = "i3.1", 
  fixed = TRUE)

# i4 <- "2017-02-14"
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-14", 
  replacement = "i3.2", 
  fixed = TRUE)

# i4 <- "2017-02-16" "2017-02-17" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-16", 
  replacement = "i4.1", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-17", 
  replacement = "i4.1", 
  fixed = TRUE)

# i6 <- "2017-02-21" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-21", 
  replacement = "i4.2", 
  fixed = TRUE)

# i7 <- "2017-02-24" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-24", 
  replacement = "i5", 
  fixed = TRUE)

# i8 <- "2017-02-28" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-28", 
  replacement = "i5", 
  fixed = TRUE)

# i9 <- "2017-03-03" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-03-03", 
  replacement = "i6", 
  fixed = TRUE)

# i10 <- "2017-03-06" "2017-03-07" "2017-03-08" "2017-03-09" "2017-03-10"
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-03-06", 
  replacement = "i6", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-03-07", 
  replacement = "i6", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-03-08", 
  replacement = "i6", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-03-09", 
  replacement = "i6", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-03-10", 
  replacement = "i6",  
  fixed = TRUE)

NROW(df)
df <- na.omit(df, cols = c("Cohort","collection_date"))
NROW(df)

head(df)

df1 <- df %>%
  select(unrooted_pd,bwpd,pc1,pc2,pc3,pc4,pc5,Cohort,collection_date,isolation_source,
         BIRTH_DAY,breed,LINE,stig,nurse)
NROW(df1)

# aggregating by avg (unique samples kept)
cols <- 1:7
df1 <- setDT(df1)[, lapply(.SD, mean), by=c(names(df1)[8:15]), .SDcols=cols]
NROW(df1)


df_breed <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$breed)$p.value,
      bwpd=kruskal.test(.$bwpd, .$breed)$p.value,
      pc1=kruskal.test(.$pc1, .$breed)$p.value,
      pc2=kruskal.test(.$pc2, .$breed)$p.value,
      pc3=kruskal.test(.$pc3, .$breed)$p.value,
      pc4=kruskal.test(.$pc4, .$breed)$p.value,
      pc5=kruskal.test(.$pc5, .$breed)$p.value,
      grouping=paste0("breed"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_breed

df_breed_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$breed)$p.value,
      bwpd=kruskal.test(.$bwpd, .$breed)$p.value,
      pc1=kruskal.test(.$pc1, .$breed)$p.value,
      pc2=kruskal.test(.$pc2, .$breed)$p.value,
      pc3=kruskal.test(.$pc3, .$breed)$p.value,
      pc4=kruskal.test(.$pc4, .$breed)$p.value,
      pc5=kruskal.test(.$pc5, .$breed)$p.value,
      grouping=paste0("breed"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_breed_all

df_line <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$LINE)$p.value,
      bwpd=kruskal.test(.$bwpd, .$LINE)$p.value,
      pc1=kruskal.test(.$pc1, .$LINE)$p.value,
      pc2=kruskal.test(.$pc2, .$LINE)$p.value,
      pc3=kruskal.test(.$pc3, .$LINE)$p.value,
      pc4=kruskal.test(.$pc4, .$LINE)$p.value,
      pc5=kruskal.test(.$pc5, .$LINE)$p.value,
      grouping=paste0("line"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_line

df_line_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$LINE)$p.value,
      bwpd=kruskal.test(.$bwpd, .$LINE)$p.value,
      pc1=kruskal.test(.$pc1, .$LINE)$p.value,
      pc2=kruskal.test(.$pc2, .$LINE)$p.value,
      pc3=kruskal.test(.$pc3, .$LINE)$p.value,
      pc4=kruskal.test(.$pc4, .$LINE)$p.value,
      pc5=kruskal.test(.$pc5, .$LINE)$p.value,
      grouping=paste0("line"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_line_all

df_bday <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_bday

# breeds
"Landrace x Cross bred (LW x D)"
"Duroc x Landrace"
"Duroc x Large white"
"Large white x Duroc"

df_bday_DurocxLandrace <- df1[df1$breed=="Duroc x Landrace",] %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day - Duroc x Landrace"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_bday_DurocxLandrace

df_bday_DurocxLw <- df1[df1$breed=="Duroc x Large white",] %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day - Duroc x Large white"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_bday_DurocxLw

df_bday_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_bday_all

df_stig <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$stig)$p.value,
      bwpd=kruskal.test(.$bwpd, .$stig)$p.value,
      pc1=kruskal.test(.$pc1, .$stig)$p.value,
      pc2=kruskal.test(.$pc2, .$stig)$p.value,
      pc3=kruskal.test(.$pc3, .$stig)$p.value,
      pc4=kruskal.test(.$pc4, .$stig)$p.value,
      pc5=kruskal.test(.$pc5, .$stig)$p.value,
      grouping=paste0("stig mother"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_stig

df_stig_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$stig)$p.value,
      bwpd=kruskal.test(.$bwpd, .$stig)$p.value,
      pc1=kruskal.test(.$pc1, .$stig)$p.value,
      pc2=kruskal.test(.$pc2, .$stig)$p.value,
      pc3=kruskal.test(.$pc3, .$stig)$p.value,
      pc4=kruskal.test(.$pc4, .$stig)$p.value,
      pc5=kruskal.test(.$pc5, .$stig)$p.value,
      grouping=paste0("stig mother"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_stig_all

df_nurse <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$nurse)$p.value,
      bwpd=kruskal.test(.$bwpd, .$nurse)$p.value,
      pc1=kruskal.test(.$pc1, .$nurse)$p.value,
      pc2=kruskal.test(.$pc2, .$nurse)$p.value,
      pc3=kruskal.test(.$pc3, .$nurse)$p.value,
      pc4=kruskal.test(.$pc4, .$nurse)$p.value,
      pc5=kruskal.test(.$pc5, .$nurse)$p.value,
      grouping=paste0("nurse mother"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_nurse

df_nurse_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$nurse)$p.value,
      bwpd=kruskal.test(.$bwpd, .$nurse)$p.value,
      pc1=kruskal.test(.$pc1, .$nurse)$p.value,
      pc2=kruskal.test(.$pc2, .$nurse)$p.value,
      pc3=kruskal.test(.$pc3, .$nurse)$p.value,
      pc4=kruskal.test(.$pc4, .$nurse)$p.value,
      pc5=kruskal.test(.$pc5, .$nurse)$p.value,
      grouping=paste0("nurse mother"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_nurse_all

df_Cohort <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$Cohort)$p.value,
      bwpd=kruskal.test(.$bwpd, .$Cohort)$p.value,
      pc1=kruskal.test(.$pc1, .$Cohort)$p.value,
      pc2=kruskal.test(.$pc2, .$Cohort)$p.value,
      pc3=kruskal.test(.$pc3, .$Cohort)$p.value,
      pc4=kruskal.test(.$pc4, .$Cohort)$p.value,
      pc5=kruskal.test(.$pc5, .$Cohort)$p.value,
      grouping=paste0("cohorts"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_Cohort

df_Cohort_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$Cohort)$p.value,
      bwpd=kruskal.test(.$bwpd, .$Cohort)$p.value,
      pc1=kruskal.test(.$pc1, .$Cohort)$p.value,
      pc2=kruskal.test(.$pc2, .$Cohort)$p.value,
      pc3=kruskal.test(.$pc3, .$Cohort)$p.value,
      pc4=kruskal.test(.$pc4, .$Cohort)$p.value,
      pc5=kruskal.test(.$pc5, .$Cohort)$p.value,
      grouping=paste0("cohorts"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_Cohort_all

all_pvalues <- rbind(df_breed_all, df_breed,
                     df_line_all, df_line, 
                     df_bday_all, df_bday, df_bday_DurocxLandrace, df_bday_DurocxLw, 
                     df_stig_all, df_stig, 
                     df_nurse_all, df_nurse, 
                     df_Cohort_all)

all_pvalues$test <- "Kruskal-Wallis"

# write out as csv
fwrite(x = all_pvalues, file = "out/all_pvalues.csv")

# write out in workbook
addWorksheet(wb, "all_pvalues")
writeData(wb, sheet = "all_pvalues", all_pvalues, rowNames = FALSE)
saveWorkbook(wb, "out/stats.xlsx", overwrite=TRUE)

padj_function <- function(x, na.rm = FALSE) (p.adjust(x,method="hommel"))

df_breed <- df_breed %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_line <- df_line %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday <- df_bday %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday_DurocxLandrace <- df_bday_DurocxLandrace %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday_DurocxLw <- df_bday_DurocxLw %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_nurse <- df_nurse %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_stig <- df_stig %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 


all_padj_Hommel <- rbind(df_breed,
                         df_line, 
                         df_bday, 
                         df_bday_DurocxLandrace, 
                         df_bday_DurocxLw, 
                         df_stig, 
                         df_nurse)


all_padj_Hommel$test <- "Kruskal-Wallis"
all_padj_Hommel$padj_method <- "Hommel"

# write out in workbook
addWorksheet(wb, "all_padj")
writeData(wb, sheet = "all_padj_Hommel", all_padj_Hommel, rowNames = FALSE)
saveWorkbook(wb, "out/stats.xlsx", overwrite=TRUE)


# adjusted pvalues

# chosen method is Tukey: 

# When you do Tukeys test, the variance is estimated from the whole set of data 
# as a pooled estimate. If the population variances are the 
# same in all groups, such a pooled estimate is much more robust and precise 
# than the individual estimated from just a part of the whole set of data. 
# Further, Tukeys procedure adjusts the p-values for multiple testing, so that 
# the family-wise error rate is controlled (probability to get at least one false 
# positive among the family of tests performed).


# by breed

aov.out = aov(unrooted_pd ~ breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$breed)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$breed)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$breed)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$breed)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$breed)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$breed)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$breed)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_breed <- rbind(aov.out1,
      aov.out2,
      aov.out3,
      aov.out4,
      aov.out5,
      aov.out6,
      aov.out7)
by_breed$group = "breed"


# by LINE

# to character otherwise considered numeric
df1$LINE <- as.character(df1$LINE)

aov.out = aov(unrooted_pd ~ LINE, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_LINE <- rbind(aov.out1,
                 aov.out2,
                 aov.out3,
                 aov.out4,
                 aov.out5,
                 aov.out6,
                 aov.out7)
by_LINE$group = "LINE"



# by BIRTH_DAY

# to character otherwise considered numeric
df1$BIRTH_DAY <- as.character(df1$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY <- rbind(aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_BIRTH_DAY$group = "BIRTH_DAY"


# by BIRTH_DAY for breed "Duroc x Landrace"

df1_sub <- df1[df1$breed=="Duroc x Landrace",]

# to character otherwise considered numeric
df1_sub$BIRTH_DAY <- as.character(df1_sub$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY_Duroc_x_Landrace <- rbind(aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_BIRTH_DAY_Duroc_x_Landrace$group = "BIRTH_DAY_Duroc_x_Landrace"


# by BIRTH_DAY for breed "Duroc x Large white"

df1_sub <- df1[df1$breed=="Duroc x Large white",]

# to character otherwise considered numeric
df1_sub$BIRTH_DAY <- as.character(df1_sub$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY_Duroc_x_Large_white <- rbind(aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_BIRTH_DAY_Duroc_x_Large_white$group = "BIRTH_DAY_Duroc_x_Large_white"


# by BIRTH_DAY for breed "Large white x Duroc"

df1_sub <- df1[df1$breed=="Large white x Duroc",]

# to character otherwise considered numeric
df1_sub$BIRTH_DAY <- as.character(df1_sub$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY_Large_white_x_Duroc <- rbind( 
                      aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_BIRTH_DAY_Large_white_x_Duroc$group = "BIRTH_DAY_Large_white_x_Duroc"


# not enough timepoint for "Landrace x Cross bred (LW x D)"


# by nurse

# to character otherwise considered numeric
df1$nurse <- as.character(df1$nurse)


aov.out = aov(unrooted_pd ~ nurse, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ nurse, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ nurse, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ nurse, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ nurse, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ nurse, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ nurse, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_nurse <- rbind( 
  aov.out1,
  aov.out2,
  aov.out3,
  aov.out4,
  aov.out5,
  aov.out6,
  aov.out7)

by_nurse$group = "nurse"

# by stig

# to character otherwise considered numeric
df1$stig <- as.character(df1$stig)


aov.out = aov(unrooted_pd ~ stig, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$stig)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ stig, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$stig)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ stig, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$stig)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ stig, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$stig)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ stig, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$stig)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ stig, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$stig)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ stig, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$stig)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_stig <- rbind( 
  aov.out1,
  aov.out2,
  aov.out3,
  aov.out4,
  aov.out5,
  aov.out6,
  aov.out7)

by_stig$group = "stig"

##################

# by Cohort

aov.out = aov(unrooted_pd ~ Cohort, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort$group = "Cohort"


# by Cohort i1

df1_sub <- df1[df1$collection_date=="i1",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i1 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i1$group = "Cohort_i1"


# by Cohort i2.1

df1_sub <- df1[df1$collection_date=="i2.1",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i2.1 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i2.1$group = "Cohort_i2.1"



# by Cohort i2.2

df1_sub <- df1[df1$collection_date=="i2.2",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i2.2 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i2.2$group = "Cohort_i2.2"



# by Cohort i3.1

df1_sub <- df1[df1$collection_date=="i3.1",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i3.1 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i3.1$group = "Cohort_i3.1"



# by Cohort i3.2

df1_sub <- df1[df1$collection_date=="i3.2",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i3.2 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i3.2$group = "Cohort_i3.2"



# by Cohort i4.1

df1_sub <- df1[df1$collection_date=="i4.1",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i4.1 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i4.1$group = "Cohort_i4.1"



# by Cohort i4.2

df1_sub <- df1[df1$collection_date=="i4.2",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i4.2 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i4.2$group = "Cohort_i4.2"



# by Cohort i5

df1_sub <- df1[df1$collection_date=="i5",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i5 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i5$group = "Cohort_i5"



# by Cohort i6

df1_sub <- df1[df1$collection_date=="i6",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_i6 <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort_i6$group = "Cohort_i6"

all_padj_Tukey <- rbind(by_breed,
             by_LINE, 
             by_BIRTH_DAY, 
             by_BIRTH_DAY_Duroc_x_Landrace,
             by_BIRTH_DAY_Duroc_x_Large_white,
             by_BIRTH_DAY_Large_white_x_Duroc, 
             by_nurse,
             by_stig,
             by_Cohort, 
             by_Cohort_i1,
             by_Cohort_i2.1,
             by_Cohort_i2.2,
             by_Cohort_i3.1,
             by_Cohort_i3.2,
             by_Cohort_i4.1,
             by_Cohort_i4.2,
             by_Cohort_i5,
             by_Cohort_i6
)

all_padj_Tukey$test <- "anova"
all_padj_Tukey$padj_method <- "TukeyHSD"

# write out in workbook
addWorksheet(wb, "all_padj_Tukey")
writeData(wb, sheet = "all_padj_Tukey", all_padj_Tukey, rowNames = FALSE)
saveWorkbook(wb, "out/stats.xlsx", overwrite=TRUE)


# plot p-values for start factors

piglets_factors <- all_pvalues %>%
  filter(grouping != "cohorts" &
           grouping != "ctrl_neo" &
           grouping != "Dscour_ColiGuard" &
           grouping != "NeoD_NeoC" &
           collection_date != "all") 

piglets_factors$grouping <- gsub("birth day - Duroc x Large white",
                                 "bday - Duroc x Large white",piglets_factors$grouping)
piglets_factors$grouping <- gsub("birth dday - Duroc x Landrace",
                                 "bday - Duroc x Landrace",piglets_factors$grouping)

piglets_factors2 <- all_padj_Hommel %>%
  filter(grouping != "cohorts" &
           grouping != "ctrl_neo" &
           grouping != "Dscour_ColiGuard" &
           grouping != "NeoD_NeoC" &
           collection_date != "all") 

piglets_factors2$grouping <- gsub("birth day - Duroc x Large white",
                                 "bday - Duroc x Large white",piglets_factors2$grouping)
piglets_factors2$grouping <- gsub("birth dday - Duroc x Landrace",
                                 "bday - Duroc x Landrace",piglets_factors2$grouping)


unroo <- ggplot(piglets_factors, aes(collection_date,unrooted_pd, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="unrooted PD - p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_point(data = piglets_factors2, color = "red", shape = 2)
bwpd <- ggplot(piglets_factors, aes(collection_date,bwpd, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="bwpd - p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_point(data = piglets_factors2, color = "red", shape = 2)
pc1 <- ggplot(piglets_factors, aes(collection_date,pc1, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc1 - p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_point(data = piglets_factors2, color = "red", shape = 2)
pc2 <- ggplot(piglets_factors, aes(collection_date,pc2, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc2 - p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_point(data = piglets_factors2, color = "red", shape = 2)
pc3 <- ggplot(piglets_factors, aes(collection_date,pc3, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc3 - p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_point(data = piglets_factors2, color = "red", shape = 2)
pc4 <- ggplot(piglets_factors, aes(collection_date,pc4, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc4 - p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_point(data = piglets_factors2, color = "red", shape = 2)
pc5 <- ggplot(piglets_factors, aes(collection_date,pc5, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc5 - p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_point(data = piglets_factors2, color = "red", shape = 2)

pdf("out/start_factors_pvalues.pdf")
ggarrange(
  unroo, bwpd,ncol=1,nrow=2, labels = c("A","B")
)
ggarrange(
  pc1, pc2,ncol=1,nrow=2, labels = c("C","D")
)
ggarrange(
  pc3, pc4,ncol=1,nrow=2, labels = c("E","F")
)
ggarrange(
  pc5,ncol=1,nrow=2, labels = c("G","H")
)
dev.off()

######################################################################################################

# rationale: maybe cohort clusters don't show up in beta div because the first 
# principal components (5 we looked at) are explaining the changes with time 
# AIM: see if clusters appear in beta diversity when running guppy on specific time intervals

# 1 # prepares metadata files to feed guppy_epca.sh
# 1 # metadata file = 1 timepoint in trial 

# 2 # run guppy_epca.sh on HPC --> output

# 3 # outputs: HPC -> local (/Users/12705859/Desktop/metapigs_base/phylosift/input_files/)
# then to R 

# 4 # remove batch effect

# 5 # merge with metadata and average out duplicate samples -> boggo 

# 6 # plot first 4 principal components 

# filter out pos controls, neg controls and mothers from metadata
mdat_sel <- mdat %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Mothers")
NROW(mdat_sel)
unique(mdat_sel$Cohort)

mdat$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"*")

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

all <- as.character(mdat_sel$ID)
a <- as.character(mdat_Ja31$ID)
b <- as.character(mdat_Fe7$ID)
c <- as.character(mdat_Fe14$ID)
d <- as.character(mdat_Fe21$ID)
e <- as.character(mdat_Fe28$ID)

writeLines(unlist(all), "mdat_all.txt", sep = " ")
writeLines(unlist(a), "mdat_Ja31.txt", sep = " ")
writeLines(unlist(b), "mdat_Fe7.txt", sep = " ")
writeLines(unlist(c), "mdat_Fe14.txt", sep = " ")
writeLines(unlist(d), "mdat_Fe21.txt", sep = " ")
writeLines(unlist(e), "mdat_Fe28.txt", sep = " ")

#settings for plots
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.title.x=element_text(colour="black",size=8),
             axis.title.y=element_text(colour="black",size=8),
             axis.text.x=element_text(colour="black",size=8),
             axis.text.y=element_text(colour="black",size=8),
             axis.ticks=element_line(colour="black"),
             legend.position="bottom",
             legend.text=element_text(size=8),
             legend.title=element_text(size=8),
             plot.margin=unit(c(0.3,0.3,0.3,0.3),"line"))

###########################################################################################

# run phylosift guppy_epca.sh

###########################################################################################

# upload pca output file

pcadat<-read_csv(paste0(basedir,"out_Ja31/pigpca.proj"),col_names = FALSE)

pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("pc1","pc2","pc3","pc4","pc5")
pcadat <- pcadat %>%
  select(DNA_plate,DNA_well,pc1,pc2,pc3,pc4,pc5)

###########################################################################################

# remove batch effect

DNA_plate <- pcadat$DNA_plate
DNA_well <- pcadat$DNA_well

pcadat<- data.matrix(pcadat[,3:7], rownames.force = NA)
pcadat_clean<-ComBat(dat=t(as.matrix(pcadat)),DNA_plate,mod=NULL)

pcadat_clean <- t(pcadat_clean)

pcadat_clean <- as.data.frame(pcadat_clean)

pcadat_clean <- unfactor(pcadat_clean[])

pcadat <- cbind(DNA_plate,DNA_well,pcadat_clean)

# checked if any batch effect still present and it seems fine
aov.out_check = aov(pc1 ~ DNA_plate, data=pcadat)   # checked for other components too
res <- TukeyHSD(aov.out_check)
aov.out_check <- as.data.frame(res$DNA_plate)
aov.out_check$`p adj`

# pcadat is now ready to use 

###########################################################################################

# 5 # merge with metadata and average out duplicate samples -> boggo 

mdat <- mdat %>%
  select(sample_name, isolation_source, collection_date, Cohort, DNA_plate, DNA_well)

# Merge alpha and div to metadata one at a time

# merge metadata with alpha div 
coggo <- inner_join(pcadat,mdat)
NROW(coggo)
unique(coggo$Cohort)

# aggregating dups 
# select the necessary columns
coggo <- coggo %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source,
         collection_date, Cohort)

# for beta we don't need to filter out then rbind the pos controls as we don't have them

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
coggo <- setDT(coggo)[, lapply(.SD, mean), by=c(names(coggo)[6:8]), .SDcols=cols]
coggo$Cohort <- factor(coggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"
                       ))

# 6 # plot first 4 principal components 

Ja31_pc1pc2 <- ggplot(coggo, aes(x=pc1,y=pc2,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)

Ja31_pc3pc4 <- ggplot(coggo, aes(x=pc3,y=pc4,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)


###########################################################################################

# upload pca output file

pcadat<-read_csv(paste0(basedir,"out_Fe7/pigpca.proj"),col_names = FALSE)

pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("pc1","pc2","pc3","pc4","pc5")
pcadat <- pcadat %>%
  select(DNA_plate,DNA_well,pc1,pc2,pc3,pc4,pc5)

###########################################################################################

# remove batch effect

DNA_plate <- pcadat$DNA_plate
DNA_well <- pcadat$DNA_well

pcadat<- data.matrix(pcadat[,3:7], rownames.force = NA)
pcadat_clean<-ComBat(dat=t(as.matrix(pcadat)),DNA_plate,mod=NULL)

pcadat_clean <- t(pcadat_clean)

pcadat_clean <- as.data.frame(pcadat_clean)

pcadat_clean <- unfactor(pcadat_clean[])

pcadat <- cbind(DNA_plate,DNA_well,pcadat_clean)

# checked if any batch effect still present and it seems fine
aov.out_check = aov(pc1 ~ DNA_plate, data=pcadat)   # checked for other components too
res <- TukeyHSD(aov.out_check)
aov.out_check <- as.data.frame(res$DNA_plate)
aov.out_check$`p adj`

# pcadat is now ready to use 

###########################################################################################

# 5 # merge with metadata and average out duplicate samples -> boggo 

mdat <- mdat %>%
  select(sample_name, isolation_source, collection_date, Cohort, DNA_plate, DNA_well)

# Merge alpha and div to metadata one at a time

# merge metadata with alpha div 
coggo <- inner_join(pcadat,mdat)
NROW(coggo)
unique(coggo$Cohort)

# aggregating dups 
# select the necessary columns
coggo <- coggo %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source,
         collection_date, Cohort)

# for beta we don't need to filter out then rbind the pos controls as we don't have them

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
coggo <- setDT(coggo)[, lapply(.SD, mean), by=c(names(coggo)[6:8]), .SDcols=cols]
coggo$Cohort <- factor(coggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"
                       ))

# 6 # plot first 4 principal components 

Fe7_pc1pc2 <- ggplot(coggo, aes(x=pc1,y=pc2,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)

Fe7_pc3pc4 <- ggplot(coggo, aes(x=pc3,y=pc4,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)


###########################################################################################

# upload pca output file

pcadat<-read_csv(paste0(basedir,"out_Fe14/pigpca.proj"),col_names = FALSE)

pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("pc1","pc2","pc3","pc4","pc5")
pcadat <- pcadat %>%
  select(DNA_plate,DNA_well,pc1,pc2,pc3,pc4,pc5)

###########################################################################################

# remove batch effect

DNA_plate <- pcadat$DNA_plate
DNA_well <- pcadat$DNA_well

pcadat<- data.matrix(pcadat[,3:7], rownames.force = NA)
pcadat_clean<-ComBat(dat=t(as.matrix(pcadat)),DNA_plate,mod=NULL)

pcadat_clean <- t(pcadat_clean)

pcadat_clean <- as.data.frame(pcadat_clean)

pcadat_clean <- unfactor(pcadat_clean[])

pcadat <- cbind(DNA_plate,DNA_well,pcadat_clean)

# checked if any batch effect still present and it seems fine
aov.out_check = aov(pc1 ~ DNA_plate, data=pcadat)   # checked for other components too
res <- TukeyHSD(aov.out_check)
aov.out_check <- as.data.frame(res$DNA_plate)
aov.out_check$`p adj`

# pcadat is now ready to use 

###########################################################################################

# 5 # merge with metadata and average out duplicate samples -> boggo 

mdat <- mdat %>%
  select(sample_name, isolation_source, collection_date, Cohort, DNA_plate, DNA_well)

# Merge alpha and div to metadata one at a time

# merge metadata with alpha div 
coggo <- inner_join(pcadat,mdat)
NROW(coggo)
unique(coggo$Cohort)

# aggregating dups 
# select the necessary columns
coggo <- coggo %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source,
         collection_date, Cohort)

# for beta we don't need to filter out then rbind the pos controls as we don't have them

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
coggo <- setDT(coggo)[, lapply(.SD, mean), by=c(names(coggo)[6:8]), .SDcols=cols]
coggo$Cohort <- factor(coggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"
                       ))

# 6 # plot first 4 principal components 

Fe14_pc1pc2 <- ggplot(coggo, aes(x=pc1,y=pc2,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)

Fe14_pc3pc4 <- ggplot(coggo, aes(x=pc3,y=pc4,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)


###########################################################################################

# upload pca output file

pcadat<-read_csv(paste0(basedir,"out_Fe21/pigpca.proj"),col_names = FALSE)

pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("pc1","pc2","pc3","pc4","pc5")
pcadat <- pcadat %>%
  select(DNA_plate,DNA_well,pc1,pc2,pc3,pc4,pc5)

###########################################################################################

# remove batch effect

DNA_plate <- pcadat$DNA_plate
DNA_well <- pcadat$DNA_well

pcadat<- data.matrix(pcadat[,3:7], rownames.force = NA)
pcadat_clean<-ComBat(dat=t(as.matrix(pcadat)),DNA_plate,mod=NULL)

pcadat_clean <- t(pcadat_clean)

pcadat_clean <- as.data.frame(pcadat_clean)

pcadat_clean <- unfactor(pcadat_clean[])

pcadat <- cbind(DNA_plate,DNA_well,pcadat_clean)

# checked if any batch effect still present and it seems fine
aov.out_check = aov(pc1 ~ DNA_plate, data=pcadat)   # checked for other components too
res <- TukeyHSD(aov.out_check)
aov.out_check <- as.data.frame(res$DNA_plate)
aov.out_check$`p adj`

# pcadat is now ready to use 

###########################################################################################

# 5 # merge with metadata and average out duplicate samples -> boggo 

mdat <- mdat %>%
  select(sample_name, isolation_source, collection_date, Cohort, DNA_plate, DNA_well)

# Merge alpha and div to metadata one at a time

# merge metadata with alpha div 
coggo <- inner_join(pcadat,mdat)
NROW(coggo)
unique(coggo$Cohort)

# aggregating dups 
# select the necessary columns
coggo <- coggo %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source,
         collection_date, Cohort)

# for beta we don't need to filter out then rbind the pos controls as we don't have them

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
coggo <- setDT(coggo)[, lapply(.SD, mean), by=c(names(coggo)[6:8]), .SDcols=cols]
coggo$Cohort <- factor(coggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"
                       ))

# 6 # plot first 4 principal components 

Fe21_pc1pc2 <- ggplot(coggo, aes(x=pc1,y=pc2,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)

Fe21_pc3pc4 <- ggplot(coggo, aes(x=pc3,y=pc4,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)


###########################################################################################

# upload pca output file

pcadat<-read_csv(paste0(basedir,"out_Fe28/pigpca.proj"),col_names = FALSE)

pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("pc1","pc2","pc3","pc4","pc5")
pcadat <- pcadat %>%
  select(DNA_plate,DNA_well,pc1,pc2,pc3,pc4,pc5)

###########################################################################################

# remove batch effect

DNA_plate <- pcadat$DNA_plate
DNA_well <- pcadat$DNA_well

pcadat<- data.matrix(pcadat[,3:7], rownames.force = NA)
pcadat_clean<-ComBat(dat=t(as.matrix(pcadat)),DNA_plate,mod=NULL)

pcadat_clean <- t(pcadat_clean)

pcadat_clean <- as.data.frame(pcadat_clean)

pcadat_clean <- unfactor(pcadat_clean[])

pcadat <- cbind(DNA_plate,DNA_well,pcadat_clean)

# checked if any batch effect still present and it seems fine
aov.out_check = aov(pc1 ~ DNA_plate, data=pcadat)   # checked for other components too
res <- TukeyHSD(aov.out_check)
aov.out_check <- as.data.frame(res$DNA_plate)
aov.out_check$`p adj`

# pcadat is now ready to use 

###########################################################################################

# 5 # merge with metadata and average out duplicate samples -> boggo 

mdat <- mdat %>%
  select(sample_name, isolation_source, collection_date, Cohort, DNA_plate, DNA_well)

# Merge alpha and div to metadata one at a time

# merge metadata with alpha div 
coggo <- inner_join(pcadat,mdat)
NROW(coggo)
unique(coggo$Cohort)

# aggregating dups 
# select the necessary columns
coggo <- coggo %>%
  select(pc1,pc2,pc3,pc4,pc5,isolation_source,
         collection_date, Cohort)

# for beta we don't need to filter out then rbind the pos controls as we don't have them

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
coggo <- setDT(coggo)[, lapply(.SD, mean), by=c(names(coggo)[6:8]), .SDcols=cols]
coggo$Cohort <- factor(coggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"
                       ))

# 6 # plot first 4 principal components 

Fe28_pc1pc2 <- ggplot(coggo, aes(x=pc1,y=pc2,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)

Fe28_pc3pc4 <- ggplot(coggo, aes(x=pc3,y=pc4,color=Cohort))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)

###########################################################################################

Ja31_pc1pc2
Ja31_pc3pc4

Fe7_pc1pc2 # D-scour separates slightly from other cohorts in pc1
Fe7_pc3pc4

Fe14_pc1pc2 # maybe ...? # Neo+D separate out slightly in pc2 # Neo+C separate out slightly in pc1 
Fe14_pc3pc4 # maybe ...? # Neo+C separate out slightly in pc4 # D-scour separate out slightly in pc4

Fe21_pc1pc2 # maybe ...? * # Neo+D separates out in pc2 from the Neomycin cohort # D-scour separate out slightly in pc2
Fe21_pc3pc4 # Neomycin smallest ellipse (80% confidence) and it sepaartes slightly from others cohorts in pc4

Fe28_pc1pc2
Fe28_pc3pc4 # maybe ...? * # Neo+C separates out in pc3 and pc4 from Neomycin
# ColiGuard separates slightly in pc4 from Control 
# D-scour separates slightly in pc3 from Control 

# saving plots that show clustering by cohort
pdf("out/cohorts_beta.pdf")
ggarrange(Fe7_pc1pc2, 
          Fe14_pc1pc2, 
          Fe14_pc3pc4,
          Fe21_pc1pc2, 
          Fe21_pc3pc4,
          Fe28_pc3pc4,
          labels=c("A","B","C","D","E","F"),
          common.legend = TRUE, legend="bottom")
dev.off()


###########################################################################################


