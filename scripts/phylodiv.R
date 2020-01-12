
#  R version 3.6.1 (2019-07-05)

# 0   # loading input data
# 1   # batch effect
# 1.1 # merge alpha div data with metadata -> boggo
# 1.2 # merge beta div data with metadata -> coggo
# 2   # plot samples distribution across plates (time and cohorts)
# 2.1 # boggo and coggo: formatting & averaging duplicates
# 3   # plot ALPHA diversity (all timepoints)
# 4   # plot BETA diversity (all timepoints)
# 5   # plot ALPHA diversity (at pig trial start) 
# 6   # plot BETA diversity (at pig trial start) 
# 7   # plot ALPHA & BETA diversity comparing cohorts during time
# 8   # p-values

###########################################################################################

# install packages

pkgs <- c("ggbiplot","ggpubr","sva","tidyverse","broom","cowplot","data.table","dunn.test","plyr",
           "dplyr","forcats","ggplot2","gridExtra","plotrix","readr","readxl","tidyr","varhandle","tibble","purr","remotes")

install.packages(pkgs[], repos='https://cran.rstudio.com', dependencies = TRUE)  

remotes::install_github("vqv/ggbiplot")

###########################################################################################

> BiocManager::install("sva")
Error: Bioconductor version '3.8' requires R version '3.5'; see
https://bioconductor.org/install



# load libraries

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

# from local 
setwd("/Users/12705859/Desktop/metapigs_base/phylosift")
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"

# from HPC
#setwd("/shared/homes/s1/pig_microbiome/metapigs_base")
#basedir = "/shared/homes/s1/pig_microbiome/metapigs_base/input_files/"

###########################################################################################

# 0   # loading input data

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

pdf("out/edge_pca_batch_all.pdf")
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

# Batch stats - before batch removal
capture.output(
  paste0("############################# batch stats - before batch removal - alpha ########### "),
  df_DNA_plate_all <- fpddat %>%
    do({
      data.frame(
        collection_date=paste0("all"),
        sample_size=NROW(.),
        unrooted_pd=kruskal.test(.$unrooted_pd, .$DNA_plate)$p.value,
        bwpd=kruskal.test(.$bwpd, .$DNA_plate)$p.value,
        grouping=paste0("DNA plate"),
        stringsAsFactors=FALSE)
    }) %>%
    ungroup(),
  df_DNA_plate_all,
  file = "out/batch_stats.txt")

capture.output(
  paste0("############################# batch stats - before batch removal - beta ########### "),
  df_DNA_plate_all <- pcadat %>%
    do({
      data.frame(
        collection_date=paste0("all"),
        sample_size=NROW(.),
        pc1=kruskal.test(.$pc1, .$DNA_plate)$p.value,
        pc2=kruskal.test(.$pc2, .$DNA_plate)$p.value,
        pc3=kruskal.test(.$pc3, .$DNA_plate)$p.value,
        pc4=kruskal.test(.$pc4, .$DNA_plate)$p.value,
        pc5=kruskal.test(.$pc5, .$DNA_plate)$p.value,
        grouping=paste0("DNA plate"),
        stringsAsFactors=FALSE)
    }) %>%
    ungroup(),
  df_DNA_plate_all,
  file = "out/batch_stats.txt",
  append = TRUE)

# shit! we have batch effects! let's look at it

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

pdf("out/alpha_beta_batch.pdf")
annotate_figure(figure,
                top = text_grob("Batch effect by alpha and beta diversity", color = "black", size = 14)
)
dev.off()


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

capture.output(
  paste0("############################# batch stats - after batch removal - alpha ########### "),
  df_DNA_plate_all <- fpddat %>%
    do({
      data.frame(
        collection_date=paste0("all"),
        sample_size=NROW(.),
        unrooted_pd=kruskal.test(.$unrooted_pd, .$DNA_plate)$p.value,
        bwpd=kruskal.test(.$bwpd, .$DNA_plate)$p.value,
        grouping=paste0("DNA plate"),
        stringsAsFactors=FALSE)
    }) %>%
    ungroup(),
  df_DNA_plate_all,
  file = "out/batch_stats.txt",
  append = TRUE)

capture.output(
  paste0("############################# batch stats - after batch removal - beta ########### "),
  df_DNA_plate_all <- pcadat %>%
    do({
      data.frame(
        collection_date=paste0("all"),
        sample_size=NROW(.),
        pc1=kruskal.test(.$pc1, .$DNA_plate)$p.value,
        pc2=kruskal.test(.$pc2, .$DNA_plate)$p.value,
        pc3=kruskal.test(.$pc3, .$DNA_plate)$p.value,
        pc4=kruskal.test(.$pc4, .$DNA_plate)$p.value,
        pc5=kruskal.test(.$pc5, .$DNA_plate)$p.value,
        grouping=paste0("DNA plate"),
        stringsAsFactors=FALSE)
    }) %>%
    ungroup(),
  df_DNA_plate_all,
  file = "out/batch_stats.txt",
  append = TRUE)


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

pdf("out/alpha_beta_no_batch.pdf")
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

pdf("out/samples_distribution_DNA_plate_timeintervals.pdf")
ggarrange(
  p1,p2,nrow=2, labels=c("A","B")
)
dev.off()

# frequency of date by Cohort

df1 <- setDT(boggo)[, .(Freq = .N), by = .(collection_date,Cohort)]

df1[order(df1$collection_date)]

pdf("out/samples_distribution_cohorts_timeintervals.pdf")
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

# 3   # plot ALPHA diversity (all timepoints)

# ALPHA diversity overall (includes pos controls):


pdf("out/alpha_phyloentropy.pdf",width=9,height=5)
par(mar=(c(5, 10, 4, 2) +0.1))
boxplot(boggo$phylo_entropy~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin+D-scour","Neomycin","Mothers","PosControl_ColiGuard","PosControl_D-scour","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Phylogenetic entropy",ylab=NULL,las=1)
dev.off()

boggo<-inner_join(fpddat,mdat)
pdf("out/alpha_unrooted.pdf",width=9,height=5)
par(mar=(c(5, 10, 4, 2) +0.1))
boxplot(boggo$unrooted_pd~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin+D-scour","Neomycin","Mothers","PosControl_ColiGuard","PosControl_D-scour","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Unrooted PD",ylab=NULL,las=1)
dev.off()

boggo<-inner_join(fpddat,mdat)
pdf("out/alpha_bwpd.pdf",width=9,height=5)
par(mar=(c(5, 10, 4, 2) +0.1))
boxplot(boggo$bwpd~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin+D-scour","Neomycin","Mothers","PosControl_ColiGuard","PosControl_D-scour","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Balance-weighted PD",ylab=NULL,las=1)
dev.off()


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

capture.output(
  paste0("################### Piglets and mothers - stats ###################"),
  summary(boggo$bwpd[boggo$Cohort!="MockCommunity"&boggo$Cohort!="NegativeControl"&boggo$Cohort!="PosControl_D-scour"&boggo$Cohort!="PosControl_ColiGuard"]),
  paste0("################### Per cohort - BWPD ###################"),
  tapply(boggo$bwpd, boggo$Cohort, summary),
  paste0("################### Per cohort - unrooted ###################"),
  tapply(boggo$unrooted_pd, boggo$Cohort, summary),
  paste0("################### Mothers - BWPD ###################"),
  summary(boggo$bwpd[boggo$Cohort=="Mothers"]),
  paste0("################### Piglets - BWPD ###################"),
  summary(boggo$bwpd[boggo$Cohort!="MockCommunity"&boggo$Cohort!="NegativeControl"&boggo$Cohort!="PosControl_D-scour"&boggo$Cohort!="PosControl_ColiGuard"&boggo$Cohort!="Mothers"]),
  paste0("################### Mothers - unrooted ###################"),
  summary(boggo$unrooted_pd[boggo$Cohort=="Mothers"]),
  paste0("################### Piglets - unrooted ###################"),
  summary(boggo$unrooted_pd[boggo$Cohort!="MockCommunity"&boggo$Cohort!="NegativeControl"&boggo$Cohort!="PosControl_D-scour"&boggo$Cohort!="PosControl_ColiGuard"&boggo$Cohort!="Mothers"]),
  paste0("################### Per cohort - mean and sd ###################"),
  ddply(boggo, "Cohort", summarise, mean=mean(unrooted_pd), sd=sd(unrooted_pd)),
  ddply(boggo, "Cohort", summarise, mean=mean(bwpd), sd=sd(bwpd)),
  file = "out/pd_numbers.txt"
)

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
pdf("out/alpha_unrooted_time.pdf",width=9,height=5)
gen_unrooted <- ggplot(doggo, aes(x=collection_date, y=unrooted_pd, 
                                  fill=Cohort)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_boxplot() + 
  labs(x = "collection date") +
  ylim(60,260) +
  stat_compare_means(comparisons = my_comparisons)+
  theme(legend.position="none")
gen_unrooted
dev.off()

# general time change - BWPD
pdf("out/alpha_bwpd_time.pdf",width=9,height=5)
gen_bwpd <- ggplot(doggo, aes(x=collection_date, y=bwpd, 
                              fill=Cohort)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_boxplot() + 
  labs(x = "collection date") +
  ylim(1.7,3.5) +
  stat_compare_means(comparisons = my_comparisons)+
  theme(legend.position="top")
gen_bwpd
dev.off()

pdf("out/alpha_unrooted&bwpd_time.pdf")
grid.arrange(
  gen_bwpd, gen_unrooted, nrow = 2
)
dev.off()

# unrooted pd - fill: collection date
pdf("out/alpha_unrooted_cohorts.pdf",width=9,height=5)
unrooted_time <- ggplot(doggo, aes(x=Cohort, y=unrooted_pd, 
                                   fill=collection_date)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot() +
  scale_fill_discrete(name = "collection date") + 
  scale_y_continuous(limits = c(60, 170))
unrooted_time
dev.off()

# bwpd - fill: collection date
pdf("out/alpha_bwpd_fill_cohorts.pdf",width=9,height=5)
bwpd_time <- ggplot(doggo, aes(x=Cohort, y=bwpd, 
                               fill=collection_date)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot() +
  scale_fill_discrete(name = "collection date") + 
  scale_y_continuous(limits = c(1.75, 2.8))
bwpd_time
dev.off()

pdf("out/alpha_unrooted&bwpd_cohorts.pdf")
grid.arrange(
  unrooted_time, bwpd_time, nrow = 2,
  top = "Alpha diversity"
)
dev.off()

#########################################################################

my_comparisons <- list( c("2017-01-31", "2017-02-07"), 
                        c("2017-02-07", "2017-02-14"), 
                        c("2017-02-07", "2017-02-21") )

pdf("out/alpha_unrooted_cohorts_facets.pdf",width=9,height=5)
p <- ggboxplot(doggo, x = "collection_date", y = "unrooted_pd",
               color = "collection_date", palette = "jco",
               add = "jitter",facet.by = "Cohort", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(0,250)
p + stat_compare_means(comparisons = my_comparisons)  # Add pairwise comparisons p-value
dev.off()
p

#########

pdf("out/alpha_bwpd_cohorts_facets.pdf",width=9,height=5)
p <- ggboxplot(doggo, x = "collection_date", y = "bwpd",
               color = "collection_date", palette = "jco",
               add = "jitter",
               facet.by = "Cohort", short.panel.labs = FALSE)+
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(1.5,3.5)
p + stat_compare_means(comparisons = my_comparisons)  # Add pairwise comparisons p-value
dev.off()


######################################################################################################

# 4   # plot BETA diversity (all timepoints)

# BETA diversity overall (does not include pos controls):

# plot

pdf("out/edge_pca.pdf")
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

rbow <- rainbow(60, end=0.7, alpha=0.7)

pdf("out/PCA_time.pdf")
plot(coggo$pc1[coggo$Cohort!="Mothers"&coggo$Cohort!="NegativeControl"],coggo$pc2[coggo$Cohort!="Mothers"&coggo$Cohort!="NegativeControl"],main="beta diversity (phylosift edge PCA)",xlab="PC1",ylab="PC2",type="p",col=rbow[as.Date(coggo$collection_date[coggo$Cohort!="Mothers"&coggo$Cohort!="NegativeControl"])-as.Date("2017-01-29 00:00:00")])
legvec <- c(0,15,30,45,60)
color_legend( -2.9, 4, 3.5, 1.5, "trial days:", legvec, rbow)
dev.off()

#################################

# timeseries within cohort 

pdf("out/PCA_cohorts.pdf")
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

pdf("out/alpha_piglets_bynurse.pdf")
ggarrange(a, b, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()


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


pdf("out/alpha_piglets_bystig.pdf")
ggarrange(c, d, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()


##################

# same but putting nurses and stigs on the same plot, dividing BWPD from unrooted
pdf("out/alpha_BWPD_bystig_bynurse.pdf")
ggarrange(c, a, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

pdf("out/alpha_unrooted_bystig_bynurse.pdf")
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
  stat_compare_means(method = "kruskal.test", label.y=80)  # Add pairwise comparisons p-value

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

pdf("out/breed.pdf")
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

pdf("out/bday.pdf")
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
  stat_compare_means(method = "kruskal.test", label.y=80)  # Add pairwise comparisons p-value

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

pdf("out/line.pdf")
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
pdf("out/unrooted_bday_bybreed.pdf",width=9,height=5)
p <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "unrooted_pd",
               color = "BIRTH_DAY", palette = "jco",
               add = "jitter",
               facet.by = "breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(50,230)
p + stat_compare_means(comparisons = my_comparisons)
dev.off()
pdf("out/bwpd_bday_bybreed.pdf",width=9,height=5)
p <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "bwpd",
               color = "BIRTH_DAY", palette = "jco",
               add = "jitter",
               facet.by = "breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(1.6,3.4)
p + stat_compare_means(comparisons = my_comparisons)
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

moms.pca1 <- prcomp(startDF11[,1:5], center = TRUE, scale. = TRUE)
moms.pca2 <- prcomp(startDF12[,1:5], center = TRUE, scale. = TRUE)
moms.pca3 <- prcomp(startDF13[,1:5], center = TRUE, scale. = TRUE)
moms.pca4 <- prcomp(startDF14[,1:5], center = TRUE, scale. = TRUE)
moms.pca5 <- prcomp(startDF15[,1:5], center = TRUE, scale. = TRUE)
moms.pca6 <- prcomp(startDF16[,1:5], center = TRUE, scale. = TRUE)

pdf("out/piglets_to_nurse1.pdf")
nurse1 <- ggbiplot(moms.pca1, obs.scale = 1, var.scale = 1,
                   groups = startDF11$nurse, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
nurse1
dev.off()
pdf("out/piglets_to_nurse2.pdf")
nurse2 <- ggbiplot(moms.pca2, obs.scale = 1, var.scale = 1,
                   groups = startDF12$nurse, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
nurse2
dev.off()
pdf("out/piglets_to_nurse3.pdf")
nurse3 <- ggbiplot(moms.pca3, obs.scale = 1, var.scale = 1,
                   groups = startDF13$nurse, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
nurse3
dev.off()
pdf("out/piglets_to_nurse4.pdf")
nurse4 <- ggbiplot(moms.pca4, obs.scale = 1, var.scale = 1,
                   groups = startDF14$nurse, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
nurse4
dev.off()
pdf("out/piglets_to_nurse5.pdf")
nurse5 <- ggbiplot(moms.pca5, obs.scale = 1, var.scale = 1,
                   groups = startDF15$nurse, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
nurse5
dev.off()
pdf("out/piglets_to_nurse6.pdf")
nurse6 <- ggbiplot(moms.pca6, obs.scale = 1, var.scale = 1,
                   groups = startDF16$nurse, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
nurse6
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

moms.pca1 <- prcomp(startDF11[,1:5], center = TRUE, scale. = TRUE)
moms.pca2 <- prcomp(startDF12[,1:5], center = TRUE, scale. = TRUE)
moms.pca3 <- prcomp(startDF13[,1:5], center = TRUE, scale. = TRUE)
moms.pca4 <- prcomp(startDF14[,1:5], center = TRUE, scale. = TRUE)
moms.pca5 <- prcomp(startDF15[,1:5], center = TRUE, scale. = TRUE)
moms.pca6 <- prcomp(startDF16[,1:5], center = TRUE, scale. = TRUE)


pdf("out/piglets_to_stig1.pdf")
stig1 <- ggbiplot(moms.pca1, obs.scale = 1, var.scale = 1,
                  groups = startDF11$stig, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
stig1
dev.off()
pdf("out/piglets_to_stig2.pdf")
stig2 <- ggbiplot(moms.pca2, obs.scale = 1, var.scale = 1,
                  groups = startDF12$stig, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
stig2
dev.off()
pdf("out/piglets_to_stig3.pdf")
stig3 <- ggbiplot(moms.pca3, obs.scale = 1, var.scale = 1,
                  groups = startDF13$stig, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
stig3
dev.off()
pdf("out/piglets_to_stig4.pdf")
stig4 <- ggbiplot(moms.pca4, obs.scale = 1, var.scale = 1,
                  groups = startDF14$stig, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
stig4
dev.off()
pdf("out/piglets_to_stig5.pdf")
stig5 <- ggbiplot(moms.pca5, obs.scale = 1, var.scale = 1,
                  groups = startDF15$stig, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
stig5
dev.off()
pdf("out/piglets_to_stig6.pdf")
stig6 <- ggbiplot(moms.pca6, obs.scale = 1, var.scale = 1,
                  groups = startDF16$stig, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'top',
        panel.background = element_rect(fill = "white", colour = "grey50"))
stig6
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

pdf("out/PCA_bars_breed.pdf")
grid.arrange(
  breed_PC1_plot, breed_PC2_plot, breed_PC3_plot, breed_PC4_plot, breed_PC5_plot, nrow = 3, ncol=2
)
dev.off()

# PCA breed plots
breed_PCA <- prcomp(startDF2[,2:6], center = TRUE, scale. = TRUE)
breedPCA12 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                       groups = startDF2$breed, ellipse = TRUE, choices=c(1,2), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
breedPCA12
breedPCA34 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                       groups = startDF2$breed, ellipse = TRUE, choices=c(3,4), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
breedPCA34
breedPCA45 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                       groups = startDF2$breed, ellipse = TRUE, choices=c(4,5), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
breedPCA45
pdf("out/PCA_dots_breed.pdf")
grid.arrange(
  breedPCA12, breedPCA34, breedPCA45, nrow = 3, ncol = 1
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

pdf("out/PCA_bars_bday.pdf")
grid.arrange(
  bday_PC1_plot, bday_PC2_plot, bday_PC3_plot, bday_PC4_plot, bday_PC5_plot, nrow = 3, ncol=2
)
dev.off()

# PCA bday plots
bday_PCA <- prcomp(startDF2[,2:6], center = TRUE, scale. = TRUE)
bdayPCA12 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                      groups = startDF2$BIRTH_DAY, ellipse = FALSE, choices=c(1,2), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
bdayPCA12
bdayPCA34 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                      groups = startDF2$BIRTH_DAY, ellipse = FALSE, choices=c(3,4), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
bdayPCA34
bdayPCA25 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                      groups = startDF2$BIRTH_DAY, ellipse = FALSE, choices=c(2,5), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
bdayPCA25

pdf("out/PCA_dots_bday.pdf")
grid.arrange(
  bdayPCA12, bdayPCA34, bdayPCA25, nrow = 3, ncol = 1
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

pdf("out/PCA_bars_line.pdf")
grid.arrange(
  line_PC1_plot, line_PC2_plot, line_PC3_plot, line_PC4_plot, line_PC5_plot, nrow = 3, ncol=2
)
dev.off()


# PCA line plots
line_PCA <- prcomp(startDF2[,2:6], center = TRUE, scale. = TRUE)
linePCA12 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                      groups = startDF2$LINE, ellipse = TRUE, choices=c(1,2), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
linePCA12
linePCA34 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                      groups = startDF2$LINE, ellipse = TRUE, choices=c(3,4), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
linePCA34
linePCA45 <- ggbiplot(breed_PCA, obs.scale = 1, var.scale = 1,
                      groups = startDF2$LINE, ellipse = TRUE, choices=c(4,5), circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = "white", colour = "grey50"))
linePCA45

pdf("out/PCA_dots_line.pdf")
grid.arrange(
  linePCA12, linePCA34, linePCA45, nrow = 3, ncol = 1
)
dev.off()


##################

# bday AND breed:

my_comparisons <- list(  c("2017-01-07", "2017-01-09"),  
                         c("2017-01-08", "2017-01-10"), 
                         c("2017-01-09", "2017-01-11"), 
                         c("2017-01-10", "2017-01-11") )

p1 <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc1",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(0,4)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'none')+
  stat_compare_means(comparisons = my_comparisons)
p2 <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc2",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(0,7)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'none')+
  stat_compare_means(comparisons = my_comparisons)
p3 <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc3",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(-2,3)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'none')+
  stat_compare_means(comparisons = my_comparisons)
p4 <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc4",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(-1,2.2)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'none')+
  stat_compare_means(comparisons = my_comparisons)
p5 <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc5",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "breed", short.panel.labs = FALSE) +
  ylim(4.7,6.6)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'none')+
  stat_compare_means(comparisons = my_comparisons)

pdf("out/PCA_bars_bdaybreeds.pdf")
grid.arrange(
  p1,p2, nrow = 1, ncol = 2
)
plot.new()
grid.arrange(
  p3,p4,p5, nrow = 1, ncol = 3
)
dev.off()


##############################################

# 7   # plot ALPHA & BETA diversity comparing cohorts during time


# join alpha and beta div DFs
# remember we have not done this before because PCA dataframe is missing some samples 
# the pos controls: 

NROW(boggo)
unique(boggo$Cohort)
NROW(coggo)
unique(coggo$Cohort)
# and 12 piglet samples

# we don't care as what follows is a comparison of the pig cohorts
finalDF <- inner_join(boggo,coggo)
NROW(finalDF)

finalDF$PigPen <- NULL

# rename the cohorts to shorter names

finalDF[11] <- lapply(
  finalDF[11], 
  gsub, 
  pattern = "Neomycin+D-scour", 
  replacement = "Neo+D", 
  fixed = TRUE)
finalDF[11] <- lapply(
  finalDF[11], 
  gsub, 
  pattern = "Neomycin+ColiGuard", 
  replacement = "Neo+C", 
  fixed = TRUE)
finalDF[11] <- lapply(
  finalDF[11], 
  gsub, 
  pattern = "Neomycin", 
  replacement = "Neo", 
  fixed = TRUE)

# reorder
finalDF$Cohort <- factor(finalDF$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neo",
                                "Neo+D",
                                "Neo+C"))

# Plot
# finalDF cohorts - finalDF values - by time interval 

my_comparisons = list( c("Control", "D-scour"), 
                        c("Control", "ColiGuard"),
                        c("Neo", "Neo+D"),
                        c("Neo", "Neo+C"),
                        c("Control", "Neo") )


#font size for pvalues 
your_font_size <- 4
# to plot the pdfs rather than the tiffs, decrease size by 2

My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 9),
  axis.text.y = element_text(size = 9),
  axis.title.y = element_text(size = 11))
# to plot the pdfs rather than the tiffs, decrease size by 2


############################### t1 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-01-31" |
    collection_date == "2017-02-01" ) 

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,250)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1.5,3.5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,4.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,6.5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,3.4)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.5,2.6)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,7.2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 4, ncol = 2
)

pdf("out/alpha_beta_t1.pdf")
annotate_figure(figure,
                top = text_grob("Trial day 2 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
ggsave("out/t1.tiff", all_plots, width=15.8, height=11)

############################### t2.1 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-02-03" ) 

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,230)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1.7,4)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,4.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,7.2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,3.4)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.5,1.4)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,7.7)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 4, ncol = 2
)

pdf("out/alpha_beta_t2.1.pdf")
annotate_figure(figure,
                top = text_grob("Trial day 5 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
ggsave("out/t2.1.tiff", all_plots, width=15.8, height=11)

############################### t2.2 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-02-07" )

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,260)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1.7,3.5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1,4.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,5.6)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,3)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.5,1.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all

# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,7.1)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all


figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 4, ncol = 2
)

pdf("out/alpha_beta_t2.2.pdf")
annotate_figure(figure,
                top = text_grob("Trial day 9 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
ggsave("out/t2.2.tiff", all_plots, width=15.8, height=11)

############################### t3.1 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-02-10" )

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,250)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
unroo
# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1.7,2.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
bw
# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-3,5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -3, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc1
# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,4.7)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc2
# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,2.2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc3
# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.5,1)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc4
# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5.5,6.7)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc5

figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 4, ncol = 2
)

pdf("out/alpha_beta_t3.1.pdf")
annotate_figure(figure,
                top = text_grob("Trial day 12 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
ggsave("out/t3.1.tiff", all_plots, width=15.8, height=11)

############################### t3.2 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-02-14" )

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(50,230)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 50, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
unroo
# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1.7,2.7)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
bw
# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc1
# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0.5,5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 0.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc2
# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  #ylim(-2,1.2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc3
# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.7,1.2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc4
# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,6.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc5

figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 4, ncol = 2
)

pdf("out/alpha_beta_t3.2.pdf")
annotate_figure(figure,
                top = text_grob("Trial day 16 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
ggsave("out/t3.2.tiff", all_plots, width=15.8, height=11)

############################### t4 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-02-16" |
    collection_date == "2017-02-17" |
    collection_date == "2017-02-21" )

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(50,220)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 50, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
unroo
# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1.7,2.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
bw
# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,4.5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc1
# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  #ylim(1,5.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc2
# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  #ylim(-2,1.2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc3
# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.5,1.3)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc4
# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,7)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc5

figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 4, ncol = 2
)

pdf("out/alpha_beta_t4.pdf")
annotate_figure(figure,
                top = text_grob("Trial days 18-23 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
ggsave("out/t4.tiff", all_plots, width=15.8, height=11)

############################### t5 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-02-24" |
    collection_date == "2017-02-28" )

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(75,180)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 75, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
unroo
# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  #ylim(1.7,2.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
bw
# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,4.7)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc1
# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1,5.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc2
# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,1.1)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc3
# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.5,2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc4
# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,7)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc5

figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 4, ncol = 2
)

pdf("out/alpha_beta_t5.pdf")
annotate_figure(figure,
                top = text_grob("Trial days 26-30 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
all_plots
ggsave("out/t5.tiff", all_plots, width=15.8, height=11)

############################### t6 ############################### 

toggo <- finalDF %>% filter(
  collection_date == "2017-03-03" |
    collection_date == "2017-03-06" |
    collection_date == "2017-03-07" |
    collection_date == "2017-03-08" |
    collection_date == "2017-03-09" |
    collection_date == "2017-03-10" )

# unrooted
unroo <- ggboxplot(toggo, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(75,180)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 75, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
unroo
# bwpd
bw <- ggboxplot(toggo, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  #ylim(1.7,2.5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 1.7, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
bw
# pc1
pc1 <- ggboxplot(toggo, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,4.2)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc1
# pc2
pc2 <- ggboxplot(toggo, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(2,5.5)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc2
# pc3
pc3 <- ggboxplot(toggo, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1,0.6)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc3
# pc4
pc4 <- ggboxplot(toggo, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  #ylim(-1.2,1.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = -1.2, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc4
# pc5
pc5 <- ggboxplot(toggo, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(5,6.8)+
  My_Theme+
  geom_hline(yintercept = mean(toggo$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(comparisons = my_comparisons, size = your_font_size) +
  stat_compare_means(method = "anova", label.y = 5, size = your_font_size) +        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE, size = your_font_size)      # Pairwise comparison against all
pc5

figure <- grid.arrange(
  unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, ncol = 4
)

pdf("out/alpha_beta_t6.pdf")
annotate_figure(figure,
                top = text_grob("Trial days 23-40 - Alpha and beta diversity among cohorts", color = "black", size = 14)
)
dev.off()

all_plots <- plot_grid(NULL, NULL, NULL, NULL, unroo, bw, pc1, pc2, pc3, pc4, pc5, nrow = 3, 
                       labels = c("A", "", "","", "B", "C", "D", "E", "F", "G", "H"),
                       ncol = 4)
ggsave("out/t6.tiff", all_plots, width=15.8, height=11)


# find out plot colors to be replicated in the timeline
scales::show_col(scales::hue_pal()(6))

######################################################################################################

# same but using facets. Looking at one measurement (unrooted, bwpd, etc)
# at a time
# separating by time intervals


df <- finalDF
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
  replacement = "i4", 
  fixed = TRUE)
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-17", 
  replacement = "i4", 
  fixed = TRUE)

# i6 <- "2017-02-21" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-21", 
  replacement = "i4", 
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

df <- na.omit(df, cols = c("Cohort","collection_date"))
NROW(df)

#font size for pvalues 
your_font_size <- 3

# other fonts
My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 7),
  axis.title.y = element_text(size = 9))

# unrooted
rep(seq(130,160,length.out=15),8) # deciding y.pos for signif values
stat.test <- df %>%
  group_by(collection_date) %>%
  t_test(unrooted_pd ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  mutate(y.position=rep(seq(130,160,length.out=15),8))
# any significant? 
NROW(which(stat.test$p.adj.signif != "ns"))
unroo <- ggboxplot(df, x = "Cohort", y = "unrooted_pd", color = "Cohort", 
                   legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(40,200)+
  facet_wrap(~collection_date)+
  My_Theme+
  geom_hline(yintercept = mean(df$unrooted_pd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 40, size = your_font_size) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = your_font_size)
unroo

# bwpd
rep(seq(2.5,3,length.out=15),8) # deciding y.pos for signif values
stat.test <- df %>%
  group_by(collection_date) %>%
  t_test(bwpd ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  mutate(y.position=rep(seq(2.5,3,length.out=15),8))
# any significant? 
NROW(which(stat.test$p.adj.signif != "ns"))
bw <- ggboxplot(df, x = "Cohort", y = "bwpd", color = "Cohort", 
                legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(1,3.3)+
  facet_wrap(~collection_date)+
  My_Theme+
  geom_hline(yintercept = mean(df$bwpd), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1, size = your_font_size) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = your_font_size)
bw

# pc1
rep(seq(2.5,4,length.out=15),8) # deciding y.pos for signif values
stat.test <- df %>%
  group_by(collection_date) %>%
  t_test(pc1 ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  mutate(y.position=rep(seq(2.5,4,length.out=15),8))
# any significant? 
NROW(which(stat.test$p.adj.signif != "ns"))
pc1 <- ggboxplot(df, x = "Cohort", y = "pc1", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,4)+
  facet_wrap(~collection_date)+
  My_Theme+
  geom_hline(yintercept = mean(df$pc1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = your_font_size)
pc1

# pc2
rep(seq(3.5,4.5,length.out=15),8) # deciding y.pos for signif values
stat.test <- df %>%
  group_by(collection_date) %>%
  t_test(pc2 ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  mutate(y.position=rep(seq(3,4,length.out=15),8))
# any significant? 
NROW(which(stat.test$p.adj.signif != "ns"))
pc2 <- ggboxplot(df, x = "Cohort", y = "pc2", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(0,4.5)+
  facet_wrap(~collection_date)+
  My_Theme+
  geom_hline(yintercept = mean(df$pc2), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 0, size = your_font_size) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = your_font_size)
pc2

# pc3
rep(seq(0.8,1.2,length.out=15),8) # deciding y.pos for signif values
stat.test <- df %>%
  group_by(collection_date) %>%
  t_test(pc3 ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  mutate(y.position=rep(seq(0.8,1.2,length.out=15),8))
# any significant? 
NROW(which(stat.test$p.adj.signif != "ns"))
pc3 <- ggboxplot(df, x = "Cohort", y = "pc3", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-2,1.5)+
  facet_wrap(~collection_date)+
  My_Theme+
  geom_hline(yintercept = mean(df$pc3), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = -2, size = your_font_size) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = your_font_size)
pc3

# pc4
rep(seq(0.5,1.2,length.out=15),8) # deciding y.pos for signif values
stat.test <- df %>%
  group_by(collection_date) %>%
  t_test(pc4 ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  mutate(y.position=rep(seq(0.5,1.2,length.out=15),8))
# any significant? 
NROW(which(stat.test$p.adj.signif != "ns"))
pc4 <- ggboxplot(df, x = "Cohort", y = "pc4", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(-1.5,1.5)+
  facet_wrap(~collection_date)+
  My_Theme+
  geom_hline(yintercept = mean(df$pc4), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = -1.5, size = your_font_size) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = your_font_size)
pc4

# pc5
rep(seq(6.3,6.8,length.out=15),8) # deciding y.pos for signif values
stat.test <- df %>%
  group_by(collection_date) %>%
  t_test(pc5 ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  mutate(y.position=rep(seq(6.3,6.8,length.out=15),8))
# any significant? 
NROW(which(stat.test$p.adj.signif != "ns"))
pc5 <- ggboxplot(df, x = "Cohort", y = "pc5", color = "Cohort", 
                 legend = "none") +
  geom_jitter(aes(colour = Cohort, x = Cohort), 
              position = position_jitter(width = .2), alpha = 1, size=0.7)+
  ylim(4.5,7)+
  facet_wrap(~collection_date)+
  My_Theme+
  geom_hline(yintercept = mean(df$pc5), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 4.5, size = your_font_size) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = your_font_size)
pc5


pdf("out/cohorts_by_interval.pdf")
grid.arrange(
  unroo, bw, nrow = 2, ncol = 1
)
grid.arrange(
  pc1,pc2, nrow = 2, ncol = 1
)
grid.arrange(
  pc3,pc4, nrow = 2, ncol = 1
)
grid.arrange(
  pc5, nrow = 2, ncol = 1
)
dev.off()

######################################################################################################

# 7.1   # p-values cohorts

df1 <- df %>%
  select(unrooted_pd,bwpd,pc1,pc2,pc3,pc4,pc5,Cohort,collection_date) %>%
  pivot_longer(cols = unrooted_pd:pc5,
               values_to = "value",
               names_to = "method")

stat.test <- df1 %>% 
  group_by(Cohort,method) %>%
  t_test(value ~ collection_date) %>%
  adjust_pvalue(method="fdr") %>%
  filter(p.adj.signif != "ns")
stat.test2 <- df1 %>%
  group_by(collection_date,method) %>%
  t_test(value ~ Cohort) %>%
  adjust_pvalue(method="fdr") %>%
  filter(p.adj.signif != "ns")

# write out 
fwrite(x = stat.test, file = "out/cohorts_WithinAndBetween_fdr_pvalues.csv")
fwrite(x = stat.test2, file = "out/cohorts_WithinAndBetween_fdr_pvalues.csv",append=TRUE)

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
  replacement = "i5.1", 
  fixed = TRUE)

# i8 <- "2017-02-28" 
df[10] <- lapply(
  df[10], 
  gsub, 
  pattern = "2017-02-28", 
  replacement = "i5.2", 
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


df_breed <- df %>%
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

df_breed_all <- df %>%
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

df_line <- df %>%
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

df_line_all <- df %>%
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

df_bday <- df %>%
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

df_bday_DurocxLandrace <- df[df$breed=="Duroc x Landrace",] %>%
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

df_bday_DurocxLw <- df[df$breed=="Duroc x Large white",] %>%
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

df_bday_all <- df %>%
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

df_stig <- df %>%
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

df_stig_all <- df %>%
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

df_nurse <- df %>%
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


df_nurse_all <- df %>%
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

df_Cohort <- df %>%
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

df_Cohort_all <- df %>%
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
df

sub <- df %>% filter(Cohort == "Control" |
                       Cohort == "Neo" ) 
df_ctrl_neo <- sub %>%
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
      grouping=paste0("ctrl_neo"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_ctrl_neo

sub <- df %>% filter(Cohort == "ColiGuard" |
                       Cohort == "D-scour" ) 
df_Dscour_ColiGuard <- sub %>%
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
      grouping=paste0("Dscour_ColiGuard"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_Dscour_ColiGuard

sub <- df %>% filter(Cohort == "Neo+D" |
                       Cohort == "Neo+C" ) 
df_neoD_neoC <- sub %>%
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
      grouping=paste0("NeoD_NeoC"),
      stringsAsFactors=FALSE)
  }) %>%
  ungroup() %>%
  select(-starts_with("i"))
df_neoD_neoC



all_pvalues <- rbind(df_breed_all, df_breed,
                     df_line_all, df_line, 
                     df_bday_all, df_bday, df_bday_DurocxLandrace, df_bday_DurocxLw, 
                     df_stig_all, df_stig, 
                     df_nurse_all, df_nurse, 
                     df_Cohort_all, df_Cohort, 
                     df_ctrl_neo, df_Dscour_ColiGuard, df_neoD_neoC)

# write out the normalized counts
fwrite(x = all_pvalues, file = "out/all_pvalues.csv")

# bh correction: 

capture.output(
  paste0("#################################### BREED ##############################"),
  dunn.test(df$unrooted_pd, 
            df$breed, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df$bwpd, 
            df$breed, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### LINE ##############################"),
  dunn.test(df$unrooted_pd, 
            df$LINE, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df$bwpd, 
            df$LINE, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### BIRTH DAY ##############################"),
  dunn.test(df$unrooted_pd, 
            df$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df$bwpd, 
            df$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### BREED ### Duroc x Landrace ############"),
  dunn.test(df[df$breed=="Duroc x Landrace",]$unrooted_pd, 
            df[df$breed=="Duroc x Landrace",]$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$breed=="Duroc x Landrace",]$bwpd, 
            df[df$breed=="Duroc x Landrace",]$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### BREED ### Duroc x Large white ############"),
  dunn.test(df[df$breed=="Duroc x Large white",]$unrooted_pd, 
            df[df$breed=="Duroc x Large white",]$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$breed=="Duroc x Large white",]$bwpd, 
            df[df$breed=="Duroc x Large white",]$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### BREED ### Large white x Duroc ############"),
  dunn.test(df[df$breed=="Large white x Duroc",]$unrooted_pd, 
            df[df$breed=="Large white x Duroc",]$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$breed=="Large white x Duroc",]$bwpd, 
            df[df$breed=="Large white x Duroc",]$BIRTH_DAY, 
            method="bh",alpha=0.05,list=TRUE),
  # not enough timepoint for "Landrace x Cross bred (LW x D)"
  paste0("#################################### Cohort ####################################"),
  dunn.test(df$unrooted_pd, 
            df$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df$bwpd, 
            df$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i1 ##################################"),
  dunn.test(df[df$collection_date=="i1",]$unrooted_pd, 
            df[df$collection_date=="i1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i1",]$bwpd, 
            df[df$collection_date=="i1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i1",]$pc1, 
            df[df$collection_date=="i1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i1",]$pc2, 
            df[df$collection_date=="i1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i1",]$pc3, 
            df[df$collection_date=="i1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i1",]$pc4, 
            df[df$collection_date=="i1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i1",]$pc5, 
            df[df$collection_date=="i1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i2.1 ##################################"),
  dunn.test(df[df$collection_date=="i2.1",]$unrooted_pd, 
            df[df$collection_date=="i2.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.1",]$bwpd, 
            df[df$collection_date=="i2.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.1",]$pc1, 
            df[df$collection_date=="i2.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.1",]$pc2, 
            df[df$collection_date=="i2.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.1",]$pc3, 
            df[df$collection_date=="i2.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.1",]$pc4, 
            df[df$collection_date=="i2.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.1",]$pc5, 
            df[df$collection_date=="i2.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i2.2 ##################################"),
  dunn.test(df[df$collection_date=="i2.2",]$unrooted_pd, 
            df[df$collection_date=="i2.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.2",]$bwpd, 
            df[df$collection_date=="i2.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.2",]$pc1, 
            df[df$collection_date=="i2.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.2",]$pc2, 
            df[df$collection_date=="i2.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.2",]$pc3, 
            df[df$collection_date=="i2.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.2",]$pc4, 
            df[df$collection_date=="i2.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i2.2",]$pc5, 
            df[df$collection_date=="i2.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i3.1 ##################################"),
  dunn.test(df[df$collection_date=="i3.1",]$unrooted_pd, 
            df[df$collection_date=="i3.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.1",]$bwpd, 
            df[df$collection_date=="i3.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.1",]$pc1, 
            df[df$collection_date=="i3.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.1",]$pc2, 
            df[df$collection_date=="i3.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.1",]$pc3, 
            df[df$collection_date=="i3.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.1",]$pc4, 
            df[df$collection_date=="i3.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.1",]$pc5, 
            df[df$collection_date=="i3.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i3.2 ##################################"),
  dunn.test(df[df$collection_date=="i3.2",]$unrooted_pd, 
            df[df$collection_date=="i3.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.2",]$bwpd, 
            df[df$collection_date=="i3.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.2",]$pc1, 
            df[df$collection_date=="i3.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.2",]$pc2, 
            df[df$collection_date=="i3.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.2",]$pc3, 
            df[df$collection_date=="i3.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.2",]$pc4, 
            df[df$collection_date=="i3.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i3.2",]$pc5, 
            df[df$collection_date=="i3.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i4.1 ##################################"),
  dunn.test(df[df$collection_date=="i4.1",]$unrooted_pd, 
            df[df$collection_date=="i4.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.1",]$bwpd, 
            df[df$collection_date=="i4.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.1",]$pc1, 
            df[df$collection_date=="i4.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.1",]$pc2, 
            df[df$collection_date=="i4.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.1",]$pc3, 
            df[df$collection_date=="i4.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.1",]$pc4, 
            df[df$collection_date=="i4.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.1",]$pc5, 
            df[df$collection_date=="i4.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i4.2 ##################################"),
  dunn.test(df[df$collection_date=="i4.2",]$unrooted_pd, 
            df[df$collection_date=="i4.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.2",]$bwpd, 
            df[df$collection_date=="i4.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.2",]$pc1, 
            df[df$collection_date=="i4.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.2",]$pc2, 
            df[df$collection_date=="i4.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.2",]$pc3, 
            df[df$collection_date=="i4.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.2",]$pc4, 
            df[df$collection_date=="i4.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i4.2",]$pc5, 
            df[df$collection_date=="i4.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i5.1 ##################################"),
  dunn.test(df[df$collection_date=="i5.1",]$unrooted_pd, 
            df[df$collection_date=="i5.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.1",]$bwpd, 
            df[df$collection_date=="i5.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.1",]$pc1, 
            df[df$collection_date=="i5.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.1",]$pc2, 
            df[df$collection_date=="i5.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.1",]$pc3, 
            df[df$collection_date=="i5.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.1",]$pc4, 
            df[df$collection_date=="i5.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.1",]$pc5, 
            df[df$collection_date=="i5.1",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i5.2 ##################################"),
  dunn.test(df[df$collection_date=="i5.2",]$unrooted_pd, 
            df[df$collection_date=="i5.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.2",]$bwpd, 
            df[df$collection_date=="i5.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.2",]$pc1, 
            df[df$collection_date=="i5.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.2",]$pc2, 
            df[df$collection_date=="i5.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.2",]$pc3, 
            df[df$collection_date=="i5.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.2",]$pc4, 
            df[df$collection_date=="i5.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i5.2",]$pc5, 
            df[df$collection_date=="i5.2",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  paste0("#################################### Cohort i6 ##################################"),
  dunn.test(df[df$collection_date=="i6",]$unrooted_pd, 
            df[df$collection_date=="i6",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i6",]$bwpd, 
            df[df$collection_date=="i6",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i6",]$pc1, 
            df[df$collection_date=="i6",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i6",]$pc2, 
            df[df$collection_date=="i6",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i6",]$pc3, 
            df[df$collection_date=="i6",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i6",]$pc4, 
            df[df$collection_date=="i6",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  dunn.test(df[df$collection_date=="i6",]$pc5, 
            df[df$collection_date=="i6",]$Cohort, 
            method="bh",alpha=0.05,list=TRUE),
  file = "out/pvalues_BH_correction.txt", append = FALSE
)

# plot p-values for start factors

piglets_factors <- all_pvalues %>%
  filter(grouping != "cohorts" &
           grouping != "ctrl_neo" &
           grouping != "Dscour_ColiGuard" &
           grouping != "NeoD_NeoC" &
           collection_date != "all") 

piglets_factors$grouping <- gsub("birth day","bday",piglets_factors$grouping)
piglets_factors$grouping <- gsub("nurse mother","nurse",piglets_factors$grouping)
piglets_factors$grouping <- gsub("stig mother","stig",piglets_factors$grouping)
piglets_factors$grouping <- gsub("bday - Duroc x Large white","bday-DxLW",piglets_factors$grouping)
piglets_factors$grouping <- gsub("bday - Duroc x Landrace","bday-DxL",piglets_factors$grouping)

unroo <- ggplot(piglets_factors, aes(fill=collection_date, y=unrooted_pd, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  labs(y="unrooted PD - p-value")+
  ylim(0,0.06)+
  theme(axis.text.x=element_text(hjust=1, angle=45),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

bwpd <- ggplot(piglets_factors, aes(fill=collection_date, y=bwpd, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  labs(y="BWPD - p-value")+
  ylim(0,0.06)+
  theme(axis.text.x=element_text(hjust=1, angle=45),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

pc1 <- ggplot(piglets_factors, aes(fill=collection_date, y=pc1, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  labs(y="PC1 - p-value")+
  ylim(0,0.06)+
  theme(axis.text.x=element_text(hjust=1, angle=45),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

pc2 <- ggplot(piglets_factors, aes(fill=collection_date, y=pc2, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  labs(y="PC3 - p-value")+
  ylim(0,0.06)+
  theme(axis.text.x=element_text(hjust=1, angle=45),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

pc3 <- ggplot(piglets_factors, aes(fill=collection_date, y=pc3, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  labs(y="PC3 - p-value")+
  ylim(0,0.06)+
  theme(axis.text.x=element_text(hjust=1, angle=45),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

pc4 <- ggplot(piglets_factors, aes(fill=collection_date, y=pc4, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  labs(y="PC4 - p-value")+
  ylim(0,0.06)+
  theme(axis.text.x=element_text(hjust=1, angle=45),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

pc5 <- ggplot(piglets_factors, aes(fill=collection_date, y=pc5, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  labs(y="PC5 - p-value")+
  ylim(0,0.06)+
  theme(axis.text.x=element_text(hjust=1, angle=45),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.position="none")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

# Extract the legend. Returns a gtable
for_legend_only <- ggplot(piglets_factors, aes(fill=collection_date, y=pc5, x=grouping)) + 
  geom_point(aes(colour = collection_date))+
  theme(axis.text.x=element_text(hjust=1, angle=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),
        legend.title = element_text(),
        legend.position="right")+
  guides(colour=guide_legend(ncol=2))
leg <- get_legend(for_legend_only)


pdf("out/piglets_start_factors_pvalues.pdf")
ggarrange(
  unroo, bwpd, pc1, pc2, pc3, pc4, pc5,leg,ncol=2,nrow=4
)
dev.off()


# just plotting again two plots per page 

piglets_factors$grouping <- gsub("bday","birth day",piglets_factors$grouping)
piglets_factors$grouping <- gsub("nurse","nurse mother",piglets_factors$grouping)
piglets_factors$grouping <- gsub("stig","stig mother",piglets_factors$grouping)
piglets_factors$grouping <- gsub("-DxLW","- Duroc x Large white",piglets_factors$grouping)
piglets_factors$grouping <- gsub("-DxL","- Duroc x Landrace",piglets_factors$grouping)

unroo <- ggplot(piglets_factors, aes(collection_date,unrooted_pd, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="unrooted PD - KW p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")
bwpd <- ggplot(piglets_factors, aes(collection_date,bwpd, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="bwpd - KW p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")
pc1 <- ggplot(piglets_factors, aes(collection_date,pc1, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc1 - KW p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")
pc2 <- ggplot(piglets_factors, aes(collection_date,pc2, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc2 - KW p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")
pc3 <- ggplot(piglets_factors, aes(collection_date,pc3, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc3 - KW p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")
pc4 <- ggplot(piglets_factors, aes(collection_date,pc4, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc4 - KW p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")
pc5 <- ggplot(piglets_factors, aes(collection_date,pc5, label = collection_date)) + 
  ylim(0,0.06)+
  labs(y="pc5 - KW p-value",
       x="")+
  geom_point(stat="identity") + 
  facet_wrap(~grouping, scales="free")+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)+
  theme(axis.text.x=element_text(size=7),
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        legend.position="none")

pdf("out/piglets_start_factors_pvalues_2plotsXpage.pdf")
ggarrange(
  unroo, bwpd,ncol=1,nrow=2
)
ggarrange(
  pc1, pc2,ncol=1,nrow=2
)
ggarrange(
  pc3, pc4,ncol=1,nrow=2
)
ggarrange(
  pc5,ncol=1,nrow=2
)
dev.off()

