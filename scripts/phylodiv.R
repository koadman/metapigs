
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
# 8   # plot distribution of cross_breeds and bdays among cohorts + cohorts distr plates 
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
#library(genefilter)
#library(compareGroups)
library(splitstackshape)
library(pheatmap) # used in pos_controls_reads.R


# from local 
setwd("/Users/12705859/Desktop/metapigs_base/phylosift")
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"

# from HPC
#setwd("/shared/homes/s1/pig_microbiome/metapigs_base")
#basedir = "/shared/homes/s1/pig_microbiome/metapigs_base/input_files/"

###########################################################################################

# 0   # loading input data

# tiffs (timelines)
timeline_deltas <- image_read(paste0(basedir,"Slide10.tiff"))
timeline_deltas_unroo <- image_read(paste0(basedir,"Slide11.tiff"))
timeline_deltas_bw <- image_read(paste0(basedir,"Slide12.tiff"))
timeline_deltas_weight <- image_read(paste0(basedir,"Slide13.tiff"))
timeline_deltas_guppy <- image_read(paste0(basedir,"Slide14.tiff"))

# load metadata 
mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)

# load details (cross_breed, line, bday, Sows)
details <- read_excel(paste0(basedir, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'pig'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse_sow'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'maternal_sow'
colnames(details)[colnames(details) == '...8'] <- 'cross_breed'
details$pig <- gsub("G","",details$pig)
details$pig <- gsub("T","",details$pig)

details <- details %>%
  dplyr::select(pig,BIRTH_DAY,LINE,cross_breed,maternal_sow,nurse_sow)
unique(details$cross_breed)


# formatting metadata column names 
mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == '*sample_name'] <- 'sample_name'

mdat <- mdat %>%
  dplyr::select(isolation_source,collection_date,Cohort,DNA_plate,DNA_well,PigPen)

# load alpha pd
fpddat<-read.table(paste0(basedir,"fpdalpha_div.tsv"),header=T,stringsAsFactors=F)
# load beta div
pcadat<-read_csv(paste0(basedir,"piggies_sel.txt.proj"),col_names = FALSE)

pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("pc1","pc2","pc3","pc4","pc5")
pcadat <- pcadat %>%
  dplyr::select(DNA_plate,DNA_well,pc1,pc2,pc3,pc4,pc5)

pcadat$DNA_plate <- gsub("plate_","P",pcadat$DNA_plate)


###########################################################################################

# create workbook to add stats 

wb <- createWorkbook()

###########################################################################################


# 1   # batch effect

# Plots the batch effect (both alpha and beta div)

fpddat$sid <- NULL

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
batch_bw <- ggboxplot(fpddat, x = "DNA_plate", y = "bwpd",
                      color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 

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

batch_pc2 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc2",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=0, size = your_font_size) 

batch_pc3 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc3",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.8, size = your_font_size) 

batch_pc4 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc4",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.5, size = your_font_size) 

batch_pc5 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=4, size = your_font_size)+
  ylim(4,7)


# Extract the legend. Returns a gtable
for_legend_only <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                             color = "DNA_plate", palette = "jco")+
  guides(color=guide_legend(ncol=3)) +
  theme(legend.position=c(0.5,0.5),  
        plot.margin=unit(c(1,1,7,1),"lines")) +
  labs(fill="") 

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

# frequency of date by Cohort

df1 <- setDT(boggo)[, .(Freq = .N), by = .(collection_date,Cohort)]

df1[order(df1$collection_date)]

p1 <- ggplot(df1, aes(fill=Cohort, y=Freq, x=collection_date)) + 
  geom_bar(position="dodge", stat="identity")+
  labs(x = "sample collection date",
       y = "number of samples",
       fill = "Cohort")+
  theme_bw()+
  theme(legend.position="top",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())

# frequency of DNA_plate by date

df1 <- setDT(boggo)[, .(Freq = .N), by = .(collection_date,DNA_plate)]

df1[order(df1$collection_date)]

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

pdf("out/distribution_cohorts_DNA_plate_time.pdf")
ggarrange(
  p1,p2,nrow=2, labels=c("A","B")
)
dev.off()


######################################################################################################

# 2.1 # boggo and coggo: formatting & averaging duplicates


# Now to get the stats and plot we need to remove (by averaging) 
# samples with identical collection date and identical isolation_source
# these are technical replicates

# aggregating dups for alpha
# select the necessary columns
boggo <- boggo %>%
  dplyr::select(phylo_entropy,quadratic,unrooted_pd,rooted_pd,bwpd,isolation_source,
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
                                "Sows",
                                "MockCommunity",
                                "PosControl_D-scour",
                                "PosControl_ColiGuard"))

# aggregating dups for beta
# select the necessary columns
coggo <- coggo %>%
  dplyr::select(pc1,pc2,pc3,pc4,pc5,isolation_source,
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
                                "Sows",
                                "NegativeControl"))




######################################################################################################

# DELTAS 


# filtering out piglets that had dysentery
boggo1 <- boggo %>%
  filter(!isolation_source=="29665"|isolation_source=="29865"|isolation_source=="29702")

pigs_1 <- boggo1 %>%
  filter(collection_date == "2017-01-31"|collection_date == "2017-02-01") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
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
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_2)

hist(pigs_2$unrooted_pd,breaks=100)
hist(pigs_2$bwpd,breaks=100)

# unrooted_pd has low values indicative of spurious sample. 
# remove those rows
pigs_2 <- pigs_2 %>%
  filter(!unrooted_pd < 50) %>%
  filter(!bwpd >2.5) 

###########################

pigs_3 <- boggo1 %>%
  filter(collection_date == "2017-02-14") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_3)

hist(pigs_3$unrooted_pd)
hist(pigs_3$bwpd)

###########################

pigs_4 <- boggo1 %>%
  filter(collection_date == "2017-02-21") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
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
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_5)

hist(pigs_5$unrooted_pd)
hist(pigs_5$bwpd)

###########################

pigs_6 <- boggo1 %>%
  filter(collection_date == "2017-03-03") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_6)

hist(pigs_6$unrooted_pd)
hist(pigs_6$bwpd)

##############################################################################

# settings for plots: 

#font size for pvalues 
your_font_size <- 2 

My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

theme_4diffs = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 9)) 

##############################################################################

# Ja31 vs Fe7

df1 <- merge(pigs_1,pigs_2, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Fe7"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe7"
res2$type <- "bwpd"
A_B <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Fe7"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe7"
res2$type <- "bwpd"
A_B_adj <- rbind(res1,res2)

#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
A_B_delta <- both
A_B_delta$time_delta <- "A_B"

##############################################################################

# Fe7 vs Fe14

df1 <- merge(pigs_2,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe7_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe14"
res2$type <- "bwpd"
B_C <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe7_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe14"
res2$type <- "bwpd"
B_C_adj <- rbind(res1,res2)


#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
B_C_delta <- both
B_C_delta$time_delta <- "B_C"


##############################################################################

# Fe14 vs Fe21

df1 <- merge(pigs_3,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe14_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe21"
res2$type <- "bwpd"
C_D <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe14_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe21"
res2$type <- "bwpd"
C_D_adj <- rbind(res1,res2)

#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
C_D_delta <- both
C_D_delta$time_delta <- "C_D"

##############################################################################

# Fe21 vs Fe28

df1 <- merge(pigs_4,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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

# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe21_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe21_vs_Fe28"
res2$type <- "bwpd"
D_E <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe21_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe21_vs_Fe28"
res2$type <- "bwpd"
D_E_adj <- rbind(res1,res2)


#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
D_E_delta <- both
D_E_delta$time_delta <- "D_E"

##############################################################################

# Fe28 vs Ma3

df1 <- merge(pigs_5,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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

# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe28_vs_Ma3"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe28_vs_Ma3"
res2$type <- "bwpd"
E_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe28_vs_Ma3"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe28_vs_Ma3"
res2$type <- "bwpd"
E_F_adj <- rbind(res1,res2)


#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
E_F_delta <- both
E_F_delta$time_delta <- "E_F"

##############################################################################

# Ja31 vs Fe14

df1 <- merge(pigs_1,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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

# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe14"
res2$type <- "bwpd"
A_C <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Fe14"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe14"
res2$type <- "bwpd"
A_C_adj <- rbind(res1,res2)


#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
A_C_delta <- both
A_C_delta$time_delta <- "A_C"


##############################################################################

# Fe7 vs Fe21

df1 <- merge(pigs_2,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe7_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe21"
res2$type <- "bwpd"
B_D <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe7_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe7_vs_Fe21"
res2$type <- "bwpd"
B_D_adj <- rbind(res1,res2)


#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
B_D_delta <- both
B_D_delta$time_delta <- "B_D"

##############################################################################

# Fe14 vs Fe28

df1 <- merge(pigs_3,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe14_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe28"
res2$type <- "bwpd"
C_E <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe14_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Fe28"
res2$type <- "bwpd"
C_E_adj <- rbind(res1,res2)


#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
C_E_delta <- both
C_E_delta$time_delta <- "C_E"


##############################################################################

# Ja31 vs Fe28

df1 <- merge(pigs_4,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe21_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe21_vs_Fe28"
res2$type <- "bwpd"
D_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe21_vs_Fe28"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe21_vs_Fe28"
res2$type <- "bwpd"
D_F_adj <- rbind(res1,res2)


#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
D_F_delta <- both
D_F_delta$time_delta <- "D_F"


##############################################################################

# Ja31 vs Fe21

df1 <- merge(pigs_1,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe21"
res2$type <- "bwpd"
A_D <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Fe21"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Fe21"
res2$type <- "bwpd"
A_D_adj <- rbind(res1,res2)


##############################################################################

# Fe14 vs Ma3

df1 <- merge(pigs_3,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe14_vs_Ma3"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Ma3"
res2$type <- "bwpd"
C_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Fe14_vs_Ma3"
res1$type <- "unrooted_pd"
res2$time_delta <- "Fe14_vs_Ma3"
res2$type <- "bwpd"
C_F_adj <- rbind(res1,res2)

#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
C_F_delta <- both
C_F_delta$time_delta <- "C_F"


##############################################################################

# Ja31 vs Ma3

df1 <- merge(pigs_1,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
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


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Ma3"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Ma3"
res2$type <- "bwpd"
A_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "Ja31_vs_Ma3"
res1$type <- "unrooted_pd"
res2$time_delta <- "Ja31_vs_Ma3"
res2$type <- "bwpd"
A_F_adj <- rbind(res1,res2)

#####################

new_df <- merge(df1,details, by.x="isolation_source",by.y="pig")
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
A_F_delta <- both
A_F_delta$time_delta <- "A_F"

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

# DELTAS p values  WITHOUT p-adjustment

# gather stats of the deltas: 
deltas_p <- rbind(A_B,B_C,C_D,D_E,E_F,A_C,B_D,C_E,D_F,A_D,C_F,A_F)
# convert rownames to first column
deltas_p <- setDT(deltas_p, keep.rownames = TRUE)[]
deltas_p$test <- "pairwise t-test"

##############################################################################

# DELTAS p values WITH p-adjustment

# gather stats of the deltas: 
deltas_padj <- rbind(A_B_adj,B_C_adj,C_D_adj,D_E_adj,E_F_adj,
                      A_C_adj,B_D_adj,C_E_adj,D_F_adj,
                      A_D_adj,C_F_adj,A_F_adj)
# convert rownames to first column
deltas_padj <- setDT(deltas_padj, keep.rownames = TRUE)[]
deltas_padj$test <- "pairwise t-test"
deltas_padj$method <- "BH"

##############################################################################

colnames(deltas_p)[colnames(deltas_p) == 'V1'] <- 'pvalues'
colnames(deltas_padj)[colnames(deltas_padj) == 'V1'] <- 'pvalues_adj'

deltas_stats <- merge(deltas_p,deltas_padj, by=c("rn","time_delta","type","test"))
deltas_stats <- cSplit(deltas_stats, "rn"," ")


deltas_stats <- deltas_stats %>%
  # join the two groups for which the p-value has been computed
  mutate(comparison=paste0(rn_1,"_vs_",rn_2)) %>%
  select(test,type,time_delta,comparison,pvalues,pvalues_adj,method)

# remove useless digits in strings "comparison"
deltas_stats$comparison <- deltas_stats$comparison %<>%
  gsub('[0-9]+', '', .)  # removes digits
  
# filtering to keep only meaningful comparisons 
# to be kept: 
meaningfulcomparisons <- c("Control_vs_ColiGuard", "ColiGuard_vs_Control",
                           "Control_vs_D-scour", "D-scour_vs_Control",
                           "Control_vs_Neomycin", "Neomycin_vs_Control",
                           "Neomycin_vs_Neomycin+D-scour", "Neomycin+D-scour_vs_Neomycin",
                           "Neomycin_vs_Neomycin+ColiGuard", "Neomycin+ColiGuard_vs_Neomycin")

# eliminate useless comparisons
deltas_stats <- deltas_stats[deltas_stats$comparison %in% meaningfulcomparisons,]


# add data to workbook 
addWorksheet(wb, "deltas_stats")
writeData(wb, sheet = "deltas_stats", deltas_stats, rowNames = FALSE)





##############################################################################
##############################################################################


# DELTAS p values WITH Tukey p-adjustment

# gather stats of the deltas: 
deltas_padj_w_model <- rbind(A_B_delta,B_C_delta,C_D_delta,D_E_delta,E_F_delta,
                             A_C_delta,B_D_delta,C_E_delta,D_F_delta,
                             C_F_delta,A_F_delta)

deltas_padj_w_model$test <- "pairwise t-test"
deltas_padj_w_model$padj_method <- "Tukey"
deltas_padj_w_model$model <- "value~Cohort+cross_breed+BIRTH_DAY"
# add data to workbook 
addWorksheet(wb, "deltas_padj_w_model")
writeData(wb, sheet = "deltas_padj_w_model", deltas_padj_w_model, rowNames = FALSE)

##############################################################################


# DELTAS - plots to visualize the proportion 
# of increase-decrease of alpha diversity by cohort
# within time frame 


##############################################################################

# settings for plots: 

#font size for pvalues 
your_font_size <- 2 

My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

theme_4diffs = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 9),
  axis.ticks.x = element_blank())

##############################################################################

# Ja31 vs Fe7

df1 <- merge(pigs_1,pigs_2, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_a_b <- df1
df_a_b$interval <- "Ja31-Fe7"


##############################################################################

# Fe7 vs Fe14

df1 <- merge(pigs_2,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_b_c <- df1
df_b_c$interval <- "Fe7-Fe14"


##############################################################################

# Fe14 vs Fe21

df1 <- merge(pigs_3,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_c_d <- df1
df_c_d$interval <- "Fe14-Fe21"

##############################################################################

# Fe21 vs Fe28

df1 <- merge(pigs_4,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_d_e <- df1
df_d_e$interval <- "Fe21-Fe28"

##############################################################################

# Ja31 vs Fe14

df1 <- merge(pigs_1,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_a_c <- df1
df_a_c$interval <- "Ja31-Fe14"

##############################################################################


# Fe7 vs Fe21

df1 <- merge(pigs_2,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_b_d <- df1
df_b_d$interval <- "Fe7-Fe21"


##############################################################################

all <- rbind(df_a_b,df_b_c, df_c_d, df_d_e,
             df_a_c, df_b_d)
head(all)


all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "Ja31-Fe7", 
  replacement = "A-B", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "Fe7-Fe14", 
  replacement = "B-C", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "Fe14-Fe21", 
  replacement = "C-D", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "Fe21-Fe28", 
  replacement = "D-E", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "Ja31-Fe14", 
  replacement = "A-C", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "Fe7-Fe21", 
  replacement = "B-D", 
  fixed = TRUE)

all$interval <- factor(all$interval,levels=c("A-B",
                                             "B-C",
                                             "C-D",
                                             "D-E",
                                             "A-C",
                                             "B-D"))



################################################################

# Group: Control, D-scour, ColiGuard, Neo

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_unroo) %>%
  mutate(.r = row_number()) 

palette <- c("#F8766D","#B79F00","#00BA38","#00BFC4")

new_df = new_df %>% 
  filter(Cohort.x=="Control"|Cohort.x=="D-scour"|Cohort.x=="ColiGuard"|Cohort.x=="Neomycin")

CTRL_unroo_deltas_facets <- ggplot(data = new_df,
                                   mapping = aes(x = .r, y = diff_unroo, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-100,50)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("unrooted PD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette)

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_bw) %>%
  mutate(.r = row_number()) # Add a row number variable

new_df = new_df %>% 
  filter(Cohort.x=="Control"|Cohort.x=="D-scour"|Cohort.x=="ColiGuard"|Cohort.x=="Neomycin")

CTRL_bw_deltas_facets <- ggplot(data = new_df,
                                mapping = aes(x = .r, y = diff_bw, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-30,20)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("BWPD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette)

################################################################

# Group: Neo, Neo+D-scour, Neo+ColiGuard

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_unroo) %>%
  mutate(.r = row_number()) 

palette_neo <- c("#00BFC4","#619CFF","#F564E3")

new_df = new_df %>% 
  filter(Cohort.x=="Neomycin"|Cohort.x=="Neomycin+D-scour"|Cohort.x=="Neomycin+ColiGuard")

NEO_unroo_deltas_facets <- ggplot(data = new_df,
                                  mapping = aes(x = .r, y = diff_unroo, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-100,50)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("unrooted PD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette_neo)

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_bw) %>%
  mutate(.r = row_number()) # Add a row number variable

new_df = new_df %>% 
  filter(Cohort.x=="Neomycin"|Cohort.x=="Neomycin+D-scour"|Cohort.x=="Neomycin+ColiGuard")

NEO_bw_deltas_facets <- ggplot(data = new_df,
                               mapping = aes(x = .r, y = diff_bw, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-30,20)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("BWPD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette_neo)

################################################################

both_bw <- plot_grid(CTRL_bw_deltas_facets,
                     NEO_bw_deltas_facets,
                     nrow=2,
                     rel_heights=c(0.5,0.5))

empty_space <- plot_grid(NULL, NULL, NULL, NULL,
                         ncol=4)

all_plots <- plot_grid(empty_space,
                       both_bw,
                       nrow=2,
                       rel_heights=c(0.25,0.5))

pdf("out/cohorts_deltas_bwpd.pdf")
ggdraw() +
  draw_image(timeline_deltas_bw, x = 0.01, y = 0.13) +
  draw_plot(all_plots)
dev.off()


################################################################

both_unroo <- plot_grid(CTRL_unroo_deltas_facets,
                        NEO_unroo_deltas_facets,
                        nrow=2,
                        rel_heights=c(0.5,0.5))

all_plots <- plot_grid(empty_space,
                       both_unroo,
                       nrow=2,
                       rel_heights=c(0.25,0.5))

pdf("out/cohorts_deltas_unroo.pdf")
ggdraw() +
  draw_image(timeline_deltas_unroo, x = 0.01, y = 0.13) +
  draw_plot(all_plots)
dev.off()

################################################################

##############################################################################

# means, standard deviations, percentage incr and decr: 

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_unroo) %>%
  dplyr::mutate(.r = row_number()) 


numbers_bw <- new_df %>% 
  group_by(interval,Cohort.x) %>% 
  dplyr::summarize(mean = mean(diff_bw),
            sum = sum(diff_bw),
            sd = sd(diff_bw),
            n = n(),
            n_pos = sum(diff_bw > 0),
            perc_pos = n_pos / n,
            n_neg = sum(diff_bw < 0),
            perc_neg = n_neg / n)

numbers_unroo <- new_df %>% 
  group_by(interval,Cohort.x) %>% 
  dplyr::summarize(mean = mean(diff_unroo),
            sum = sum(diff_unroo),
            sd = sd(diff_unroo),
            n = n(),
            n_pos = sum(diff_unroo > 0),
            perc_pos = n_pos / n,
            n_neg = sum(diff_unroo < 0),
            perc_neg = n_neg / n)

numbers_all_bw <- new_df %>% 
  group_by(interval) %>% 
  dplyr::summarize(mean = mean(diff_bw),
            sum = sum(diff_bw),
            sd = sd(diff_bw),
            n = n(),
            n_pos = sum(diff_bw > 0),
            perc_pos = n_pos / n,
            n_neg = sum(diff_bw < 0),
            perc_neg = n_neg / n)

numbers_all_unroo <- new_df %>% 
  group_by(interval) %>% 
  dplyr::summarize(mean = mean(diff_unroo),
            sum = sum(diff_unroo),
            sd = sd(diff_unroo),
            n = n(),
            n_pos = sum(diff_unroo > 0),
            perc_pos = n_pos / n,
            n_neg = sum(diff_unroo < 0),
            perc_neg = n_neg / n)


numbers_unroo$type = "unrooted_pd"
numbers_bw$type = "bwpd"
numbers_all_unroo$type = "unrooted_pd"
numbers_all_bw$type = "bwpd"
numbers_all_unroo$type = "unrooted_pd"
numbers_all_unroo$Cohort.x = "all"
numbers_all_bw$Cohort.x = "all"

numbers_unroo <- as.data.frame(numbers_unroo)
numbers_all_unroo <- as.data.frame(numbers_all_unroo)
numbers_bw <- as.data.frame(numbers_bw)
numbers_all_bw <- as.data.frame(numbers_all_bw)


# into a single dataframe : 

numbers <- rbind(numbers_all_unroo,numbers_unroo,numbers_all_bw,numbers_bw)

numbers$dates <- numbers$interval

numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "Ja31-Fe7", 
  pattern = "A-B", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "Fe7-Fe14", 
  pattern = "B-C", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "Fe14-Fe21", 
  pattern = "C-D", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "Fe21-Fe28", 
  pattern = "D-E", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "Ja31-Fe14", 
  pattern = "A-C", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "Fe7-Fe21", 
  pattern = "B-D", 
  fixed = TRUE)


deltas_percent_change <- numbers 

# add data to workbook 
addWorksheet(wb, "deltas_percent_change")
writeData(wb, sheet = "deltas_percent_change", deltas_percent_change, rowNames = FALSE)


######################################################################################################

# 3   # plot ALPHA diversity (all timepoints)

# ALPHA diversity overall (includes pos controls):

boggo1 <- boggo %>%
  filter(!isolation_source == "NegativeControl")

# reordering
boggo1$Cohort <- factor(boggo1$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard",
                                "Sows",
                                "MockCommunity",
                                "PosControl_D-scour",
                                "PosControl_ColiGuard"))

# boxplots again for alpha diversity, to be plotted in the same pdf
p1 <- ggplot(boggo1, aes(x=Cohort, y=phylo_entropy)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()
p2 <- ggplot(boggo1, aes(x=Cohort, y=bwpd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()
p3 <- ggplot(boggo1, aes(x=Cohort, y=unrooted_pd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()


pdf("out/alpha_BWPD_unrooted.pdf")
plot_grid(p2,p3,nrow=2, labels = c("A","B"))
dev.off()

pdf("out/alpha_all.pdf")
plot_grid(p1,p2,p3,nrow=3)
dev.off()

##############################################

a <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Sows") %>% 
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
  filter(!Cohort=="Sows") %>% 
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


e <- boggo %>% 
  na.omit(.) %>%
  group_by(Cohort) %>% 
  dplyr::summarise(min = min(unrooted_pd)
            ,max = max(unrooted_pd)
            ,mean = mean(unrooted_pd)
            ,n = n()
            ,sd = sd(unrooted_pd)
            ,q25 = quantile(unrooted_pd, .25)
            ,q75 = quantile(unrooted_pd, .75))

f <- boggo %>% 
  na.omit(.) %>%
  group_by(Cohort) %>% 
  dplyr::summarise(min = min(bwpd)
            ,max = max(bwpd)
            ,mean = mean(bwpd)
            ,n = n()
            ,sd = sd(bwpd)
            ,q25 = quantile(bwpd, .25)
            ,q75 = quantile(bwpd, .75)) 

a$Cohort="piglets"
b$Cohort="Sows"
c$Cohort="piglets"
d$Cohort="Sows"

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

# alpha_timeseries_all_cohorts multiple time points (6)

doggo <- boggo
doggo$collection_date <- as.character(doggo$collection_date)

doggo <- doggo %>% filter(collection_date == "2017-01-31" |
                            collection_date == "2017-02-07" |
                            collection_date == "2017-02-14" |
                            collection_date == "2017-02-21" |
                            collection_date == "2017-02-28" |
                            collection_date == "2017-03-03") 

# filter out heavy outliers

doggo <- doggo %>% 
  filter(!unrooted_pd < 25)

doggo <- doggo %>% 
  filter(!bwpd > 2.7)

# check
hist(doggo$unrooted_pd)
hist(doggo$bwpd)

# (no need to filter out pos, neg controls and Sows 
# as filtering by date is already doing it )

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
                                         "2017-02-21",
                                         "2017-02-28",
                                         "2017-03-03"))

my_comparisons = list( c("2017-01-31", "2017-02-07"), 
                       c("2017-02-07", "2017-02-14"), 
                       c("2017-02-14", "2017-02-21"),
                       c("2017-02-07", "2017-02-21"),
                       c("2017-02-21", "2017-02-28"),
                       c("2017-02-14", "2017-02-28"),
                       c("2017-02-28", "2017-03-03"))

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
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=0.1,
                size=0.2,
                position=position_dodge(0.5)) +
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, size = 2)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=0.1,
                size=0.2,
                position=position_dodge(0.5)) +
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, size=2)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "BWPD - mean")
gen_bwpd


tosave <- ggarrange(gen_unrooted, gen_bwpd, ncol = 2, labels=c("A","B"), common.legend = TRUE)
ggsave(filename = "out/time_alpha.pdf", plot = tosave)


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

unroo_all <- doggo %>%
  group_by(collection_date) %>%
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

bw_all <- doggo %>%
  group_by(collection_date) %>%
  dplyr::summarise(min = min(bwpd)
                   ,max = max(bwpd)
                   ,mean = mean(bwpd)
                   ,sd = sd(bwpd)
                   ,n = n()
                   ,q25 = quantile(bwpd, .25)
                   ,q75 = quantile(bwpd, .75)) 

unroo_all$Cohort <- "all"
bw_all$Cohort <- "all"

# add data to workbook 

bw$type="bwpd"
bw_all$type="bwpd"
unroo$type="unrooted_pd"
unroo_all$type="unrooted_pd"
bw <- as.data.frame(bw)
bw_all <- as.data.frame(bw_all)
unroo <- as.data.frame(unroo)
unroo_all <- as.data.frame(unroo_all)
both <- rbind(bw,unroo,bw_all,unroo_all,means)
addWorksheet(wb, "alpha_means")
writeData(wb, sheet = "alpha_means", both, rowNames = FALSE)



######################################################################################################


# alpha diversity change during time in cohorts, with unrooted PD & BWPD in the same plot:

x <- doggo %>%
  select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd) %>%
  droplevels()

summs <- x %>% group_by(Cohort,collection_date) %>% 
  dplyr::summarise(`mean BWPD` = mean(bwpd),
                   sdUnroo = sd(unrooted_pd),
                   meanUnroo = mean(unrooted_pd)) %>%
  droplevels()

summs <- as.data.frame(summs)

pdf("out/time_alpha_cohorts_unroo_BW.pdf")
ggplot(summs, aes(x=collection_date, y=meanUnroo,group=Cohort,color=Cohort)) +
  geom_errorbar(aes(ymin=meanUnroo-sdUnroo, ymax=meanUnroo+sdUnroo), 
                width=0.1,
                size=0.2,
                position=position_dodge(0.5)) +
  geom_line() + geom_point(aes(size=`mean BWPD`)) +
  facet_wrap(~Cohort)+
  labs(y="mean unrooted PD")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))
dev.off()

#####################################################

# powerpoint "Desktop/metapgs_base/phylosift/PD_change_cohorts_time.pptx"
# 
# xx <- summs %>%
#   dplyr::mutate(rankedUnroo = rank(meanUnroo))
# 
# # number of groups large enough to appreciate diffs between timepoints 
# xx$groupUnroo <- as.numeric(cut_number(xx$rankedUnroo,4))
# 
# Relabel = c(4,6,8,10)
# xx$groupUnroo = Relabel[xx$groupUnroo]
# 
# # BWPD into 4 categories
# 
# # number of groups large enough to appreciate diffs between timepoints 
# xx$groupBW <- as.numeric(cut_number(xx$meanBW,4))
# 
# Relabel = c(0.25,0.50,0.75,1.00)
# xx$groupBW = Relabel[xx$groupBW]
# xx

######################################################################################################


# alpha_timeseries_all_cohorts 4 time points 

doggo <- boggo
doggo$collection_date <- as.character(doggo$collection_date)

doggo <- doggo %>% filter(collection_date == "2017-01-31" |
                            collection_date == "2017-02-07" |
                            collection_date == "2017-02-14" |
                            collection_date == "2017-02-21" ) 


doggo[3] <- lapply(
  doggo[3], 
  gsub, 
  pattern = "Neomycin+D-scour", 
  replacement = "Neo+D-scour", 
  fixed = TRUE)

doggo[3] <- lapply(
  doggo[3], 
  gsub, 
  pattern = "Neomycin+ColiGuard", 
  replacement = "Neo+ColiGuard", 
  fixed = TRUE)

# filter out heavy outliers

doggo <- doggo %>% 
  filter(!unrooted_pd < 25)

doggo <- doggo %>% 
  filter(!bwpd > 2.7)

# check
hist(doggo$unrooted_pd)
hist(doggo$bwpd)

# (no need to filter out pos, neg controls and Sows 
# as filtering by date is already doing it )

# reordering
doggo$Cohort <- factor(doggo$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neo+D-scour",
                                "Neo+ColiGuard"))

# reordering
doggo$collection_date <- factor(doggo$collection_date,
                                levels=c("2017-01-31",
                                         "2017-02-07",
                                         "2017-02-14",
                                         "2017-02-21"))

my_comparisons = list( c("2017-01-31", "2017-02-07"), 
                       c("2017-02-07", "2017-02-14"), 
                       c("2017-02-14", "2017-02-21"))

stat.test_unroo <- doggo %>%
  group_by(Cohort) %>%
  t_test(unrooted_pd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") %>%
  dplyr::mutate(y.position=rep(seq(150,220,length.out=6),6)) %>%
  dplyr::mutate_if(is.numeric, round, digits = 4)

unroo <- ggboxplot(doggo, x = "collection_date", y = "unrooted_pd",
                   color = "collection_date", palette = "jco",
                   add = "jitter",short.panel.labs = TRUE) +
  theme_bw()+
  facet_grid(~Cohort)+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  ylim(50,225)+
  stat_pvalue_manual(stat.test_unroo, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = 2)

#########

stat.test_bwpd <- doggo %>%
  group_by(Cohort) %>%
  t_test(bwpd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") %>%
  mutate(y.position=rep(seq(2.3,2.65,length.out=6),6)) %>%
  mutate_if(is.numeric, round, digits = 4)

bw <- ggboxplot(doggo, x = "collection_date", y = "bwpd",
                color = "collection_date", palette = "jco",
                add = "jitter",short.panel.labs = TRUE) +
  facet_grid(~Cohort)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  ylim(1.5,2.70)+
  stat_pvalue_manual(stat.test_bwpd, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = 2)
bw

tosave <- ggarrange(unroo,bw,nrow=2,ncol=1,labels=c("A","B"),common.legend = TRUE)
ggsave(filename = "out/time_alpha_cohorts.pdf", plot = tosave)


# add data to workbook 
both <- rbind(stat.test_bwpd,
              stat.test_unroo)
both$padj_method <- "bonferroni"

addWorksheet(wb, "alpha_cohorts")
writeData(wb, sheet = "alpha_cohorts", both, rowNames = FALSE)


#######################

# comparing timepoints, irrespective of cohort

stat.test_unroo_all <- doggo %>%
  t_test(unrooted_pd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") 

stat.test_bwpd_all <- doggo %>%
  t_test(bwpd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") 

# add data to workbook 
both <- rbind(stat.test_bwpd_all,
              stat.test_unroo_all)
both$padj_method <- "bonferroni"

addWorksheet(wb, "alpha_time")
writeData(wb, sheet = "alpha_time", both, rowNames = FALSE)

######################################################################################################
######################################################################################################

# 5   # plot ALPHA diversity (at pig trial start) 

# ALPHA diversity in piglets at the start of the trial 

startDF <- boggo 

startDF <- startDF %>% filter(
  collection_date == "2017-01-31" |
    collection_date == "2017-02-01")  %>%
  dplyr::select(phylo_entropy,quadratic,unrooted_pd,rooted_pd,bwpd,isolation_source)

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
  dplyr::select(phylo_entropy,unrooted_pd,bwpd,isolation_source,nurse_sow,maternal_sow,
                cross_breed,BIRTH_DAY,LINE)


# plots

cw_summary <- startDF1 %>% 
  group_by(nurse_sow) %>% 
  tally()

a <- ggplot(startDF1, aes(x=nurse_sow, y=bwpd, group=nurse_sow)) + 
  labs(title = "Piglets alpha diversity (BWPD)",
       subtitle = "Grouped by nurse sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(nurse_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11))

b <- ggplot(startDF1, aes(x=nurse_sow, y=unrooted_pd, group=nurse_sow)) + 
  labs(title = "Piglets alpha diversity (unrooted)",
       subtitle = "Grouped by nurse sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(nurse_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11)) +
  scale_y_continuous(limits=c(60,160))

#######

cw_summary <- startDF1 %>% 
  group_by(maternal_sow) %>% 
  tally()

c <- ggplot(startDF1, aes(x=maternal_sow, y=bwpd, group=maternal_sow)) + 
  labs(title = "Piglets alpha diversity (BWPD)",
       subtitle = "Grouped by maternal sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(maternal_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11))

d <- ggplot(startDF1, aes(x=maternal_sow, y=unrooted_pd, group=maternal_sow)) + 
  labs(title = "Piglets alpha diversity (unrooted PD)",
       subtitle = "Grouped by maternal sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(maternal_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11)) +
  scale_y_continuous(limits=c(60,160))

##################

# nurse_sows and maternal sows on the same plot, dividing BWPD from unrooted
pdf("out/nursesow_maternalsow_BWPD.pdf")
ggarrange(c, a, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

pdf("out/nursesow_maternalsow_unrooted.pdf")
ggarrange(d, b, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

##############################################

# do piglets at arrival cluster by unrooted and bwpd
# based on cross_breed/line/birth day? 

startDF1$BIRTH_DAY <- as.character(startDF1$BIRTH_DAY)
startDF1$BIRTH_DAY <- factor(startDF1$BIRTH_DAY, 
                             levels=c("2017-01-06", 
                                      "2017-01-07", 
                                      "2017-01-08",
                                      "2017-01-09",
                                      "2017-01-10",
                                      "2017-01-11"))
startDF1$LINE <- as.character(startDF1$LINE)

# plots

# by cross_breed

####### get sample size within each cross_breed group:

cw_summary <- startDF1 %>% 
  group_by(cross_breed) %>% 
  tally()

# cross_breed - unrooted 
cross_breed_unrooted_plot <- ggboxplot(startDF1, x = "cross_breed", y = "unrooted_pd",
                                 color = "cross_breed", palette = "jco",
                                 add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(50,160)+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# cross_breed - bwpd

cross_breed_bwpd_plot <- ggboxplot(startDF1, x = "cross_breed", y = "bwpd",
                             color = "cross_breed", palette = "jco",
                             add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(),legend.position="none")+
  ylim(1.5,2.7)+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward")+
  stat_compare_means(method = "kruskal.test", label.y=1.5) 

tosave <- ggarrange(cross_breed_unrooted_plot, cross_breed_bwpd_plot, nrow = 2, labels=c("A","B"))
ggsave(file = "out/breed_alpha.pdf", tosave)


# by birth day

####### get sample size within each bday group:

cw_summary <- startDF1 %>% 
  group_by(BIRTH_DAY) %>% 
  tally()

# bday - unrooted 
bday_unrooted_plot <- ggboxplot(startDF1, x = "BIRTH_DAY", y = "unrooted_pd",
                                color = "BIRTH_DAY", palette = "jco",
                                add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(50,160)+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# bday - bwpd 
bday_bwpd_plot <- ggboxplot(startDF1, x = "BIRTH_DAY", y = "bwpd",
                            color = "BIRTH_DAY", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  ylim(1.5,2.7)+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=1.5)  # Add pairwise comparisons p-value


tosave <- ggarrange(bday_unrooted_plot, bday_bwpd_plot, nrow = 2, labels=c("A","B"))
ggsave(file = "out/bday_alpha.pdf", tosave)


# by line

####### get sample size within each bday group:

cw_summary <- startDF1 %>% 
  group_by(LINE) %>% 
  tally()

# line - unrooted 
LINE_unrooted_plot <- ggboxplot(startDF1, x = "LINE", y = "unrooted_pd",
                                color = "LINE", palette = "jco",
                                add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(50,160)+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# line - bwpd 
LINE_bwpd_plot <- ggboxplot(startDF1, x = "LINE", y = "bwpd",
                            color = "LINE", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  ylim(1.5,2.75)+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=1.5)  # Add pairwise comparisons p-value


tosave <- ggarrange(LINE_unrooted_plot, LINE_bwpd_plot, nrow = 2, labels=c("A","B"))
ggsave(file = "out/line_alpha.pdf", tosave)


##################

# bday AND cross_breed:
# with age unrooted pd decreases and bwpd increases 
# conclusion comes from two cross_breeds as each bday is not represented 
# equally by the 4 cross_breeds 
my_comparisons <- list(  c("2017-01-07", "2017-01-09"),  
                         c("2017-01-08", "2017-01-10"), 
                         c("2017-01-09", "2017-01-11"), 
                         c("2017-01-10", "2017-01-11") )

startDF1_sub <- startDF1 %>%
  filter(!cross_breed=="Landrace x Cross bred (LW x D)")
startDF1_sub <- startDF1_sub %>%
  filter(!cross_breed=="Large white x Duroc")


p1 <- ggboxplot(startDF1_sub, x = "BIRTH_DAY", y = "unrooted_pd",
               color = "BIRTH_DAY", palette = "jco",
               add = "jitter",
               facet.by = "cross_breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(50,230)+
  stat_compare_means(comparisons = my_comparisons)

p2 <- ggboxplot(startDF1_sub, x = "BIRTH_DAY", y = "bwpd",
               color = "BIRTH_DAY", palette = "jco",
               add = "jitter",
               facet.by = "cross_breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(1.6,3.4)+
  stat_compare_means(comparisons = my_comparisons)


tosave <- ggarrange(p1,p2, nrow = 2, labels=c("A","B"))
ggsave(file = "out/bday_bybreed_alpha.pdf", tosave)


######################################################################################################

# 6   # plot BETA diversity (at pig trial start) 

# BETA diversity in piglets at the start of the trial 

startDF <- coggo %>% filter(
  collection_date == "2017-01-31" |
    collection_date == "2017-02-01")  %>%
  dplyr::select(pc1,pc2,pc3,pc4,pc5,isolation_source)

startDF1 <- merge(startDF,details, by.x="isolation_source",by.y="pig")

startDF1 <- startDF1 %>%
  dplyr::select(pc1,pc2,pc3,pc4,pc5,isolation_source,nurse_sow,maternal_sow)

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

# nurse_sows

startDF1_unique <- startDF1 %>% group_by(isolation_source) %>% slice(1)

# remove rows were nurse_sow is unique (can't be plotted a pca)
startDF1_unique <- subset(startDF1_unique,duplicated(nurse_sow) | duplicated(nurse_sow, fromLast=TRUE))

# order alphabetically by nurse_sow
startDF1_unique <- startDF1_unique[order(startDF1_unique$nurse_sow),]

rownames(startDF1_unique) <- 1:nrow(startDF1_unique)
length(unique(startDF1_unique$nurse_sow))

# subsetting to have 5 nurse_sows in each 
startDF11 <- startDF1_unique[1:27,]
startDF12 <- startDF1_unique[28:47,]
startDF13 <- startDF1_unique[48:63,]
startDF14 <- startDF1_unique[64:80,]
startDF15 <- startDF1_unique[81:98,]
startDF16 <- startDF1_unique[99:122,]

p11<-ggplot(startDF11,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p12<-ggplot(startDF12,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p13<-ggplot(startDF13,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p14<-ggplot(startDF14,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p15<-ggplot(startDF15,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p16<-ggplot(startDF16,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, 
               level = 0.80)
p16

pdf("out/nursesow_PC1PC2.pdf")
grid.arrange(p11,p12,p13,p14,p15,p16, nrow=3,ncol=2)
dev.off()

######################

# maternal sows 

startDF1_unique <- startDF1 %>% group_by(isolation_source) %>% slice(1)

# remove rows were nurse_sow is unique (can't be plotted a pca)
startDF1_unique <- subset(startDF1_unique,duplicated(maternal_sow) | duplicated(maternal_sow, fromLast=TRUE))

# order alphabetically by nurse_sow
startDF1_unique <- startDF1_unique[order(startDF1_unique$maternal_sow),]
rownames(startDF1_unique) <- 1:nrow(startDF1_unique)
length(unique(startDF1_unique$maternal_sow))

# subsetting to have 5 maternal sows in each 
startDF11 <- startDF1_unique[1:30,]
startDF12 <- startDF1_unique[31:58,]
startDF13 <- startDF1_unique[59:72,]
startDF14 <- startDF1_unique[73:85,]
startDF15 <- startDF1_unique[86:108,]
startDF16 <- startDF1_unique[109:123,]

p11<-ggplot(startDF11,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p12<-ggplot(startDF12,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p13<-ggplot(startDF13,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p14<-ggplot(startDF14,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p15<-ggplot(startDF15,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)
p16<-ggplot(startDF16,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme+
  stat_ellipse(inherit.aes = TRUE, 
               level = 0.80)
p16

pdf("out/maternalsow_PC1PC2.pdf")
grid.arrange(p11,p12,p13,p14,p15,p16, nrow=3,ncol=2)
dev.off()

##############################################

# do piglets at arrival cluster by PCA
# based on cross_breed/line/birth day? 

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


# by cross_breed

####### get sample size within each cross_breed group:

cw_summary <- startDF2 %>% 
  group_by(cross_breed) %>% 
  tally()

# cross_breed - PCA
cross_breed_PC1_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc1",
                            color = "cross_breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.5)  # Add pairwise comparisons p-value
cross_breed_PC1_plot
cross_breed_PC2_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc2",
                            color = "cross_breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=1)  # Add pairwise comparisons p-value
cross_breed_PC2_plot
cross_breed_PC3_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc3",
                            color = "cross_breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.5)  # Add pairwise comparisons p-value
cross_breed_PC3_plot
cross_breed_PC4_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc4",
                            color = "cross_breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.2)  # Add pairwise comparisons p-value
cross_breed_PC4_plot
cross_breed_PC5_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc5",
                            color = "cross_breed", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value
cross_breed_PC5_plot


#tosave <- ggarrange(cross_breed_PC1_plot, cross_breed_PC2_plot, cross_breed_PC3_plot, cross_breed_PC4_plot, 
#          cross_breed_PC5_plot, nrow = 3, ncol=2, labels = c("A","B","C","D","E"),
#          common.legend = TRUE)
#ggsave(file = "out/cross_breed_beta.pdf", tosave)

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
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.2)  # Add pairwise comparisons p-value
bday_PC1_plot
bday_PC2_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc2",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=1)  # Add pairwise comparisons p-value
bday_PC2_plot
bday_PC3_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc3",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.5)  # Add pairwise comparisons p-value
bday_PC3_plot
bday_PC4_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc4",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1)  # Add pairwise comparisons p-value
bday_PC4_plot
bday_PC5_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc5",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value
bday_PC5_plot


tosave <- ggarrange(bday_PC1_plot, bday_PC2_plot, bday_PC3_plot, bday_PC4_plot, 
                    bday_PC5_plot, nrow = 3, ncol=2, labels = c("A","B","C","D","E"),
                    common.legend = TRUE)
ggsave(file = "out/bday_beta.pdf", tosave)


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
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.2)  # Add pairwise comparisons p-value
line_PC1_plot
line_PC2_plot <- ggboxplot(startDF2, x = "LINE", y = "pc2",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=.2)  # Add pairwise comparisons p-value
line_PC2_plot
line_PC3_plot <- ggboxplot(startDF2, x = "LINE", y = "pc3",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.4)  # Add pairwise comparisons p-value
line_PC3_plot
line_PC4_plot <- ggboxplot(startDF2, x = "LINE", y = "pc4",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1)  # Add pairwise comparisons p-value
line_PC4_plot
line_PC5_plot <- ggboxplot(startDF2, x = "LINE", y = "pc5",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value
line_PC5_plot


# tosave <- ggarrange(line_PC1_plot, line_PC2_plot, line_PC3_plot, line_PC4_plot, 
#                     line_PC5_plot,nrow = 3, ncol=2, labels = c("A","B","C","D","E"),
#                     common.legend = TRUE)
# ggsave(file = "out/line_beta.pdf", tosave)


##################

# bday AND cross_breed:

my_comparisons <- list(  c("2017-01-07", "2017-01-09"),  
                         c("2017-01-08", "2017-01-10"), 
                         c("2017-01-09", "2017-01-11"), 
                         c("2017-01-10", "2017-01-11") )

# need to exclude two cross_breeds as these two cross_breeds don't have enough
# age groups to be plotted and compared 
startDF2_sub <- startDF2 %>%
  filter(!cross_breed=="Landrace x Cross bred (LW x D)")
startDF2_sub <- startDF2_sub %>%
  filter(!cross_breed=="Large white x Duroc")


p1 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc1",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(0,4)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p2 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc2",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(0,7)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p3 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc3",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(-2,3)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p4 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc4",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(-1,2.2)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p5 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc5",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
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

# distribution of cross_breeds, birth days across cohorts

finalDF <- inner_join(boggo,coggo)
NROW(finalDF)

df <- merge(finalDF,details, by.x="isolation_source",by.y="pig")
head(df)

# distribution of cross_breeds and bdays across cohorts
df1 <- setDT(df)[, .(Freq = .N), by = .(BIRTH_DAY,cross_breed,Cohort)]
df1[order(df1$cross_breed)]

p1 <- ggplot(df1, aes(fill=BIRTH_DAY, y=Freq, x=Cohort)) + 
  geom_bar(position="stack", stat="identity")+
  theme(legend.position="none")+
  facet_wrap(.~cross_breed,scales="free")
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


######################################################################################################

# 8   # p-values

# PVALUES - 


# merge alpha&beta div (finalDF) to details and details metadata

df <- merge(finalDF,details, by.x="isolation_source",by.y="pig")
head(df)
class(df$collection_date)

unique(df$collection_date)
# rename collection dates to time intervals 
# iM <- "2017-01-30"
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-01-30", 
  replacement = "iM", 
  fixed = TRUE)
# i0 <- "2017-01-31" "2017-02-01" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-01-31", 
  replacement = "i1", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-01", 
  replacement = "i1", 
  fixed = TRUE)

# i1 <- "2017-02-03" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-03", 
  replacement = "i2.1", 
  fixed = TRUE)

# i2 <- "2017-02-06" "2017-02-07" "2017-02-08"
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-06", 
  replacement = "i2.2", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-07", 
  replacement = "i2.2", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-08", 
  replacement = "i2.2", 
  fixed = TRUE)

# i3 <- "2017-02-10" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-10", 
  replacement = "i3.1", 
  fixed = TRUE)

# i4 <- "2017-02-14"
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-14", 
  replacement = "i3.2", 
  fixed = TRUE)

# i4 <- "2017-02-16" "2017-02-17" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-16", 
  replacement = "i4.1", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-17", 
  replacement = "i4.1", 
  fixed = TRUE)

# i6 <- "2017-02-21" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-21", 
  replacement = "i4.2", 
  fixed = TRUE)

# i7 <- "2017-02-24" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-24", 
  replacement = "i5", 
  fixed = TRUE)

# i8 <- "2017-02-28" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-02-28", 
  replacement = "i5", 
  fixed = TRUE)

# i9 <- "2017-03-03" 
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-03-03", 
  replacement = "i6", 
  fixed = TRUE)

# i10 <- "2017-03-06" "2017-03-07" "2017-03-08" "2017-03-09" "2017-03-10"
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-03-06", 
  replacement = "i6", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-03-07", 
  replacement = "i6", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-03-08", 
  replacement = "i6", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-03-09", 
  replacement = "i6", 
  fixed = TRUE)
df[2] <- lapply(
  df[2], 
  gsub, 
  pattern = "2017-03-10", 
  replacement = "i6",  
  fixed = TRUE)

NROW(df)
df <- na.omit(df, cols = c("Cohort","collection_date"))
NROW(df)

head(df)

df1 <- df %>%
  dplyr::select(unrooted_pd,bwpd,pc1,pc2,pc3,pc4,pc5,Cohort,collection_date,isolation_source,
         BIRTH_DAY,cross_breed,LINE,maternal_sow,nurse_sow)
NROW(df1)

# for some reasons df1$pc2 is character and not numeric. convert: 
df1$pc2 <- as.numeric(df1$pc2)

# aggregating by avg (unique samples kept)
cols <- 1:7
df1 <- setDT(df1)[, lapply(.SD, mean), by=c(names(df1)[8:15]), .SDcols=cols]
NROW(df1)

df_cross_breed <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$cross_breed)$p.value,
      bwpd=kruskal.test(.$bwpd, .$cross_breed)$p.value,
      pc1=kruskal.test(.$pc1, .$cross_breed)$p.value,
      pc2=kruskal.test(.$pc2, .$cross_breed)$p.value,
      pc3=kruskal.test(.$pc3, .$cross_breed)$p.value,
      pc4=kruskal.test(.$pc4, .$cross_breed)$p.value,
      pc5=kruskal.test(.$pc5, .$cross_breed)$p.value,
      grouping=paste0("cross_breed"),
      stringsAsFactors=FALSE)
  })
df_cross_breed

df_cross_breed_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$cross_breed)$p.value,
      bwpd=kruskal.test(.$bwpd, .$cross_breed)$p.value,
      pc1=kruskal.test(.$pc1, .$cross_breed)$p.value,
      pc2=kruskal.test(.$pc2, .$cross_breed)$p.value,
      pc3=kruskal.test(.$pc3, .$cross_breed)$p.value,
      pc4=kruskal.test(.$pc4, .$cross_breed)$p.value,
      pc5=kruskal.test(.$pc5, .$cross_breed)$p.value,
      grouping=paste0("cross_breed"),
      stringsAsFactors=FALSE)
  }) 
df_cross_breed_all

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
  })
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
  }) 
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
  })
df_bday

# cross_breeds
"Landrace x Cross bred (LW x D)"
"Duroc x Landrace"
"Duroc x Large white"
"Large white x Duroc"

df_bday_DurocxLandrace <- df1[df1$cross_breed=="Duroc x Landrace",] %>%
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
  })
df_bday_DurocxLandrace

df_bday_DurocxLw <- df1[df1$cross_breed=="Duroc x Large white",] %>%
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
  })
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
  }) 
df_bday_all

df_maternal_sow <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$maternal_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$maternal_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$maternal_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$maternal_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$maternal_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$maternal_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$maternal_sow)$p.value,
      grouping=paste0("maternal_sow"),
      stringsAsFactors=FALSE)
  }) 
df_maternal_sow

df_maternal_sow_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$maternal_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$maternal_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$maternal_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$maternal_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$maternal_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$maternal_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$maternal_sow)$p.value,
      grouping=paste0("maternal_sow"),
      stringsAsFactors=FALSE)
  }) 
df_maternal_sow_all

df_nurse_sow <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$nurse_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$nurse_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$nurse_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$nurse_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$nurse_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$nurse_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$nurse_sow)$p.value,
      grouping=paste0("nurse_sow"),
      stringsAsFactors=FALSE)
  })
df_nurse_sow

df_nurse_sow_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$nurse_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$nurse_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$nurse_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$nurse_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$nurse_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$nurse_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$nurse_sow)$p.value,
      grouping=paste0("nurse_sow"),
      stringsAsFactors=FALSE)
  })
df_nurse_sow_all

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
  })
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
  })
df_Cohort_all

df_cross_breed_all <- as.data.frame(df_cross_breed_all)
df_cross_breed <- as.data.frame(df_cross_breed)
df_line_all <- as.data.frame(df_line_all)
df_line <- as.data.frame(df_line)
df_bday_all <- as.data.frame(df_bday_all)
df_bday <- as.data.frame(df_bday)
df_bday_DurocxLandrace <- as.data.frame(df_bday_DurocxLandrace)
df_bday_DurocxLw <- as.data.frame(df_bday_DurocxLw)
df_maternal_sow_all <- as.data.frame(df_maternal_sow_all)
df_maternal_sow <- as.data.frame(df_maternal_sow)
df_nurse_sow_all <- as.data.frame(df_nurse_sow_all)
df_nurse_sow <- as.data.frame(df_nurse_sow)
df_Cohort_all <- as.data.frame(df_Cohort_all)

all_pvalues <- rbind(df_cross_breed_all, df_cross_breed,
                     df_line_all, df_line, 
                     df_bday_all, df_bday, df_bday_DurocxLandrace, df_bday_DurocxLw, 
                     df_maternal_sow_all, df_maternal_sow, 
                     df_nurse_sow_all, df_nurse_sow, 
                     df_Cohort_all)

all_pvalues$test <- "Kruskal-Wallis"


# write out in workbook
addWorksheet(wb, "all_pvalues")
writeData(wb, sheet = "all_pvalues", all_pvalues, rowNames = FALSE)

padj_function <- function(x, na.rm = FALSE) (p.adjust(x,method="hommel"))

df_cross_breed <- df_cross_breed %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_line <- df_line %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday <- df_bday %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday_DurocxLandrace <- df_bday_DurocxLandrace %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday_DurocxLw <- df_bday_DurocxLw %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_nurse_sow <- df_nurse_sow %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_maternal_sow <- df_maternal_sow %>%
  mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 


all_padj_Hommel <- rbind(df_cross_breed,
                         df_line, 
                         df_bday, 
                         df_bday_DurocxLandrace, 
                         df_bday_DurocxLw, 
                         df_maternal_sow, 
                         df_nurse_sow)


all_padj_Hommel$test <- "Kruskal-Wallis"
all_padj_Hommel$padj_method <- "Hommel"

# write out in workbook
addWorksheet(wb, "all_padj_Hommel")
writeData(wb, sheet = "all_padj_Hommel", all_padj_Hommel, rowNames = FALSE)


# adjusted pvalues

# chosen method is Tukey: 

# When you do Tukeys test, the variance is estimated from the whole set of data 
# as a pooled estimate. If the population variances are the 
# same in all groups, such a pooled estimate is much more robust and precise 
# than the individual estimated from just a part of the whole set of data. 
# Further, Tukeys procedure adjusts the p-values for multiple testing, so that 
# the family-wise error rate is controlled (probability to get at least one false 
# positive among the family of tests performed).


# by cross_breed

aov.out = aov(unrooted_pd ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_cross_breed <- rbind(aov.out1,
      aov.out2,
      aov.out3,
      aov.out4,
      aov.out5,
      aov.out6,
      aov.out7)
by_cross_breed$group = "cross_breed"


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


# by BIRTH_DAY for cross_breed "Duroc x Landrace"

df1_sub <- df1[df1$cross_breed=="Duroc x Landrace",]

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


# by BIRTH_DAY for cross_breed "Duroc x Large white"

df1_sub <- df1[df1$cross_breed=="Duroc x Large white",]

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


# by BIRTH_DAY for cross_breed "Large white x Duroc"

df1_sub <- df1[df1$cross_breed=="Large white x Duroc",]

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


# by nurse_sow

# to character otherwise considered numeric
df1$nurse_sow <- as.character(df1$nurse_sow)


aov.out = aov(unrooted_pd ~ nurse_sow, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_nurse_sow <- rbind( 
  aov.out1,
  aov.out2,
  aov.out3,
  aov.out4,
  aov.out5,
  aov.out6,
  aov.out7)

by_nurse_sow$group = "nurse_sow"

# by maternal_sow

# to character otherwise considered numeric
df1$maternal_sow <- as.character(df1$maternal_sow)


aov.out = aov(unrooted_pd ~ maternal_sow, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_maternal_sow <- rbind( 
  aov.out1,
  aov.out2,
  aov.out3,
  aov.out4,
  aov.out5,
  aov.out6,
  aov.out7)

by_maternal_sow$group = "maternal_sow"

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

all_padj_Tukey <- rbind(by_cross_breed,
             by_LINE, 
             by_BIRTH_DAY, 
             by_BIRTH_DAY_Duroc_x_Landrace,
             by_BIRTH_DAY_Duroc_x_Large_white,
             by_BIRTH_DAY_Large_white_x_Duroc, 
             by_nurse_sow,
             by_maternal_sow,
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


# plot p-values for start factors

piglets_factors <- all_pvalues %>%
  dplyr::filter(grouping != "cohorts" &
           grouping != "ctrl_neo" &
           grouping != "Dscour_ColiGuard" &
           grouping != "NeoD_NeoC" &
           collection_date != "all") 

piglets_factors$grouping <- gsub("birth day - Duroc x Large white",
                                 "birth day - DxLW",piglets_factors$grouping)
piglets_factors$grouping <- gsub("birth day - Duroc x Landrace",
                                 "birth day - DxL",piglets_factors$grouping)

piglets_factors2 <- all_padj_Hommel %>%
  filter(grouping != "cohorts" &
           grouping != "ctrl_neo" &
           grouping != "Dscour_ColiGuard" &
           grouping != "NeoD_NeoC" &
           collection_date != "all") 

piglets_factors2$grouping <- gsub("birth day - Duroc x Large white",
                                 "birth day - DxLW",piglets_factors2$grouping)
piglets_factors2$grouping <- gsub("birth day - Duroc x Landrace",
                                 "birth day - DxL",piglets_factors2$grouping)


unique(piglets_factors$grouping)
unique(piglets_factors2$grouping)

# order to show facets:
piglets_factors$grouping <- factor(piglets_factors$grouping,
                    levels=c("birth day",
                             "birth day - DxL",
                             "birth day - DxLW",
                             "cross_breed",
                             "maternal_sow",
                             "nurse_sow",
                             "line"))
piglets_factors2$grouping <- factor(piglets_factors2$grouping,
                                 levels=c("birth day",
                                          "birth day - DxL",
                                          "birth day - DxLW",
                                          "cross_breed",
                                          "maternal_sow",
                                          "nurse_sow",
                                          "line"))

# alpha
df <- piglets_factors %>%
  pivot_longer(cols=c('unrooted_pd','bwpd'))
colnames(df)[colnames(df)=="name"] <- "parameter"
df2 <- piglets_factors2 %>%
  pivot_longer(cols=c('unrooted_pd','bwpd'))
colnames(df2)[colnames(df2)=="name"] <- "parameter"
alpha_plot <- ggplot(df, aes(x=collection_date,y=value)) + 
  ylim(0,0.06)+
  labs(y="alpha diversity - p-value",
       x="")+
  geom_point(data=df2,aes(shape=parameter), color="red", size=2)+
  geom_point(aes(shape=parameter), color="black", size=2)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="right")+
  facet_wrap(~grouping)+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)
# beta 
df <- piglets_factors %>%
  pivot_longer(cols=c(contains('pc')))
colnames(df)[colnames(df)=="name"] <- "parameter"
df2 <- piglets_factors2 %>%
  pivot_longer(cols=c(contains('pc')))
colnames(df2)[colnames(df2)=="name"] <- "parameter"
beta_plot <- ggplot(df, aes(x=collection_date,y=value)) + 
  ylim(0,0.06)+
  labs(y="beta diversity - p-value",
       x="")+
  geom_point(data=df2,aes(shape=parameter), color="red", size=2)+
  geom_point(aes(shape=parameter), color="black", size=2)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="right")+
  facet_wrap(~grouping)+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)



pdf("out/start_factors_pvalues.pdf")
ggarrange(
  alpha_plot, beta_plot,ncol=1,nrow=2, labels = c("A","B")
)
dev.off()




######################################################################################################
######################################################################################################
######################################################################################################


# save stats in workbook
saveWorkbook(wb, "/Users/12705859/Desktop/metapigs_base/phylosift/out/stats.xlsx", overwrite=TRUE)


###########################################################################################

# cite packages 

sink("/Users/12705859/Desktop/metapigs_base/metapigs_base_packages_citations.bib")
out <- sapply(names(sessionInfo()$otherPkgs), 
              function(x) print(citation(x), style = "Bibtex"))

closeAllConnections()

