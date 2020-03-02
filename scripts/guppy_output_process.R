
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
library(COMBAT)

# manual settings 

removebatcheffect_allowed <- "yes"


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
  dplyr::select(isolation_source,BIRTH_DAY,breed,stig,nurse)


###########################################################################################


my.basedir <-"~/Desktop/metapigs_base/phylosift/guppy/guppy_output"

jplace_files = list.files(my.basedir,pattern=".proj")

# construct an empty dataframe to build on 
jplace_df <- data.frame(
  DNA_plate = character(),
  DNA_well = character(),
  file = character(),
  PC1 = character(),
  PC2 = character(),
  PC3 = character(),
  PC4 = character(),
  PC5 = character(),
  stringsAsFactors = FALSE
  )

for (jplace_file in jplace_files) {
  
  # read in file 
  pcadat <- read_csv(file.path(my.basedir,jplace_file), col_names = FALSE)

  pcadat <- cSplit(pcadat, "X1","_")
  
  pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
  pcadat$DNA_well <- pcadat$X1_3
  colnames(pcadat)[1:5] <- c("PC1","PC2","PC3","PC4","PC5")
  pcadat <- pcadat %>%
    dplyr::select(DNA_plate,DNA_well,PC1,PC2,PC3,PC4,PC5)
  
  pcadat$file <- basename(jplace_file)
  
  jplace_df <- rbind(
    jplace_df, 
    pcadat
  )
  
}

# convert PC columns to numeric class 
jplace_df <- jplace_df %>%
  mutate_at('PC1',as.numeric) %>% 
  mutate_at('PC2',as.numeric) %>% 
  mutate_at('PC3',as.numeric) %>% 
  mutate_at('PC4',as.numeric) %>% 
  mutate_at('PC5',as.numeric) 

# clean the file names
jplace_df$file <- gsub('_sel.txt.proj', '', jplace_df$file)


##############################

# function to adjust p-value
padj_function <- function(x, na.rm = FALSE) (p.adjust(x,method="hommel"))

# determine if and where there is a batch effect
checkbatch_before <- jplace_df %>%
  group_by(file) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      PC1=kruskal.test(.$PC1, .$DNA_plate)$p.value,
      PC2=kruskal.test(.$PC2, .$DNA_plate)$p.value,
      PC3=kruskal.test(.$PC3, .$DNA_plate)$p.value,
      PC4=kruskal.test(.$PC4, .$DNA_plate)$p.value,
      PC5=kruskal.test(.$PC5, .$DNA_plate)$p.value,
      batch_removal=paste0("before"),
      stringsAsFactors=FALSE)
  }) %>%
  mutate_at(c("PC1","PC2","PC3","PC4","PC5"),padj_function) 

# picking of groups is based on (adj) p-value of batch effect being <0.05
batch_affected <- checkbatch_before %>%
  filter(PC1<0.05|PC2<0.05|PC3<0.05|PC4<0.05|PC5<0.05) %>%
  dplyr::select(file)

to_remove_batch <- jplace_df[jplace_df$file %in% batch_affected$file,]

# splitting into multiple dataframes (by file name)
multiple_DFs <- split( to_remove_batch , f = to_remove_batch$file )

# construct an empty dataframe to build on 
unbatched <- data.frame(
    DNA_plate = character(),
    DNA_well = character(),
    file = character(),
    PC1 = character(),
    PC2 = character(),
    PC3 = character(),
    PC4 = character(),
    PC5 = character(),
    stringsAsFactors = FALSE
  )
  

##############################
# this loop is entered if manually allowed (top of script)

if (removebatcheffect_allowed=="yes") {
  for (single_DF in multiple_DFs) {
    
    DNA_plate <- single_DF$DNA_plate
    DNA_well <- single_DF$DNA_well
    file <- single_DF$file
  
    single_DF<- data.matrix(single_DF[,4:8], rownames.force = NA)
    
    single_DF<-ComBat(dat=t(as.matrix(single_DF)),DNA_plate,mod=NULL)
  
    single_DF <- t(single_DF)
    single_DF <- as.data.frame(single_DF)
    single_DF <- unfactor(single_DF[])
    single_DF <- cbind(DNA_plate,DNA_well,file,single_DF)

    unbatched <- rbind(
    unbatched, 
    single_DF
    )
  }
} else {
  print("No batch effect removal allowed")
}

##############################

# check batch effect AFTER batch effect removal 
checkbatch_unbatched <- unbatched %>%
  group_by(file) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      PC1=kruskal.test(.$PC1, .$DNA_plate)$p.value,
      PC2=kruskal.test(.$PC2, .$DNA_plate)$p.value,
      PC3=kruskal.test(.$PC3, .$DNA_plate)$p.value,
      PC4=kruskal.test(.$PC4, .$DNA_plate)$p.value,
      PC5=kruskal.test(.$PC5, .$DNA_plate)$p.value,
      batch_removal=paste0("before"),
      stringsAsFactors=FALSE)
  }) %>%
  mutate_at(c("PC1","PC2","PC3","PC4","PC5"),padj_function) 


##############################


# re-join the dataframes (original with unbatched one)
  
rest <- anti_join(jplace_df,unbatched,by=c("DNA_plate","DNA_well","file"))
jplace_df_final <- rbind(rest,unbatched)

NROW(jplace_df_final)


##############################


# Time to plot! 


###########################################################################################


#settings for plots
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             plot.title = element_text(),
             axis.title.x=element_text(colour="black",size=8),
             axis.title.y=element_text(colour="black",size=8),
             axis.text.x=element_text(colour="black",size=8),
             axis.text.y=element_text(colour="black",size=8),
             axis.ticks=element_line(colour="black"),
             legend.position="top",
             plot.margin=unit(c(0.3,0.3,0.3,0.3),"line"))

color_legend <- function(x, y, xlen, ylen, main, tiks, colors){
  text(x, y+.6, main, adj=c(0,0), cex=1.3)
  color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.8)
}
rbow <- rainbow(40, end=0.7, alpha=0.7)

##############################


# add one level of grouping (e.g.: all group_A* files, all timepoints, belong together)

jplace_df_final$group <- jplace_df_final$file
jplace_df_final$group <- gsub('piggies_group_A', 'groupA', jplace_df_final$group)
jplace_df_final$group <- gsub('piggies_group_B', 'groupB', jplace_df_final$group)
# note to self: ^ and $ for the exact (and entire) string match 
jplace_df_final$group <- gsub('^piggies$', 'piggies_all', jplace_df_final$group)
jplace_df_final <- cSplit(jplace_df_final, "group","_")
jplace_df_final$group_1 <- gsub('pos', 'positive_controls', jplace_df_final$group_1)
jplace_df_final$group_2 <- gsub('controls', 'all_replicates', jplace_df_final$group_2)
jplace_df_final <- setnames(jplace_df_final, old = c('group_1','group_2'), new = c('sample_type','guppied_date'))
head(jplace_df_final)

##############################

# merge metadata with details (breed,bday,nurse,...)
mdat_deets <- left_join(mdat,details)

# merge metadata with beta diversity data
multi_coggo <- inner_join(jplace_df_final,mdat_deets, by=c("DNA_plate","DNA_well"))
multi_coggo <- multi_coggo %>%
  dplyr::select(DNA_plate,DNA_well,sample_name,isolation_source,collection_date,Cohort,breed,BIRTH_DAY,nurse,stig,sample_type,guppied_date,PC1,PC2,PC3,PC4,PC5)


##############################

# give some order to the variables 

multi_coggo$sample_type <- factor(multi_coggo$sample_type, 
                      levels=c("positive_controls","piggies","groupA","groupB"))

multi_coggo$guppied_date <- factor(multi_coggo$guppied_date, 
                                      levels=c("all_replicates",
                                               "all",
                                               "Ja31",
                                               "Fe7",
                                               "Fe14",
                                               "Fe21",
                                               "Fe28",
                                               "Ma3"))

multi_coggo$BIRTH_DAY <- factor(multi_coggo$BIRTH_DAY, 
                       levels=c("2017-01-06", 
                                "2017-01-07", 
                                "2017-01-08",
                                "2017-01-09",
                                "2017-01-10",
                                "2017-01-11"))

##############################

# splitting into multiple dataframes (by sample_type name)

multi_coggo <- split( multi_coggo , f = multi_coggo$sample_type )


##############################


# get dataframes 


##############################
##############################

DF_positive_controls <- as.data.frame(multi_coggo$positive_controls)

##############################
##############################

DF_piggies <- as.data.frame(multi_coggo$piggies)
DF_piggies <- DF_piggies %>%
  filter(sample_type=="piggies") %>%
  filter(guppied_date=="all")

##############################
##############################

DF_piggies_time <- as.data.frame(multi_coggo$piggies)
DF_piggies_time_Ja31 <- DF_piggies_time %>%
  filter(sample_type=="piggies") %>%
  filter(guppied_date=="Ja31")
DF_piggies_time_Fe7 <- DF_piggies_time %>%
  filter(sample_type=="piggies") %>%
  filter(guppied_date=="Fe7")
DF_piggies_time_Fe14 <- DF_piggies_time %>%
  filter(sample_type=="piggies") %>%
  filter(guppied_date=="Fe14")
DF_piggies_time_Fe21 <- DF_piggies_time %>%
  filter(sample_type=="piggies") %>%
  filter(guppied_date=="Fe21")
DF_piggies_time_Fe28 <- DF_piggies_time %>%
  filter(sample_type=="piggies") %>%
  filter(guppied_date=="Fe28")
DF_piggies_time_Ma3 <- DF_piggies_time %>%
  filter(sample_type=="piggies") %>%
  filter(guppied_date=="Ma3")
DF_piggies_time <- rbind(
  DF_piggies_time_Ja31,
  DF_piggies_time_Fe7,
  DF_piggies_time_Fe14,
  DF_piggies_time_Fe21,
  DF_piggies_time_Fe28,
  DF_piggies_time_Ma3)

##############################
##############################

groupA <- as.data.frame(multi_coggo$groupA)
groupA_Ja31 <- groupA %>%
  filter(guppied_date=="Ja31")
groupA_Fe7 <- groupA %>%
  filter(guppied_date=="Fe7")
groupA_Fe14 <- groupA %>%
  filter(guppied_date=="Fe14")
groupA_Fe21 <- groupA %>%
  filter(guppied_date=="Fe21")
groupA_Fe28 <- groupA %>%
  filter(guppied_date=="Fe28")
groupA_Ma3 <- groupA %>%
  filter(guppied_date=="Ma3")
groupA <- rbind(
  groupA_Ja31,
  groupA_Fe7,
  groupA_Fe14,
  groupA_Fe21,
  groupA_Fe28,
  groupA_Ma3)

##############################
##############################

groupB <- as.data.frame(multi_coggo$groupB)
groupB_Ja31 <- groupB %>%
  filter(guppied_date=="Ja31")
groupB_Fe7 <- groupB %>%
  filter(guppied_date=="Fe7")
groupB_Fe14 <- groupB %>%
  filter(guppied_date=="Fe14")
groupB_Fe21 <- groupB %>%
  filter(guppied_date=="Fe21")
groupB_Fe28 <- groupB %>%
  filter(guppied_date=="Fe28")
groupB_Ma3 <- groupB %>%
  filter(guppied_date=="Ma3")
groupB <- rbind(
  groupB_Ja31,
  groupB_Fe7,
  groupB_Fe14,
  groupB_Fe21,
  groupB_Fe28,
  groupB_Ma3)

##############################
##############################


# PLOT! 


##############################
##############################

# positive controls

DF_positive_controls$Cohort <- factor(DF_positive_controls$Cohort, 
                       levels=c("MockCommunity",
                                "PosControl_D-scour",
                                "PosControl_ColiGuard"))

PC1PC2_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC2,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.5,0.5,4,4),"cm"))
PC3PC4_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC3,y=PC4,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE)+
  theme(plot.margin=unit(c(0.5,1.5,1.5,1.7),"cm"))
PC1PC5_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC5,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE)+
  theme(plot.margin=unit(c(0.5,0.5,1.9,1.7),"cm"))

##############################


xmldata <- pos_controls 

pdf("pos_controls.pdf")
### plot PC1PC2
PC1PC2_pos_controls
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC1(xmldata)),"%"), x = unit(0.65, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC2(xmldata)),"%"), x = unit(0.1, "npc"), 
          y = unit(0.55, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
### plot PC3PC4 
PC3PC4_pos_controls
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC3(xmldata)),"%"), x = unit(0.65, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC4(xmldata)),"%"), x = unit(0.1, "npc"), 
          y = unit(0.55, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
### plot PC1PC5 
PC1PC5_pos_controls
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC1(xmldata)),"%"), x = unit(0.65, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC5
grid.text(PC_down(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC5(xmldata)),"%"), x = unit(0.1, "npc"), 
          y = unit(0.55, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
dev.off()




##############################
##############################


# piggies (all time points)

DF_piggies # for plots
xmldata <- piggies



pdf("piggies.pdf")
# plot PC1PC2
plot(DF_piggies$PC1,DF_piggies$PC2,
     main="guppy on all - change with time",
     xlab=paste0("PC1  ",get_var(find_PC1(xmldata)),"%"),ylab=paste0("PC2  ",get_var(find_PC2(xmldata)),"%"),
     type="p",cex=0.8,
     col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.2, "npc"), 
          y = unit(0.05, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.8, "npc"), 
          y = unit(0.05, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.05, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.05, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))

# plot PC3PC4
plot(DF_piggies$PC3,DF_piggies$PC4,
     main="guppy on all - change with time",
     xlab=paste0("PC3  ",get_var(find_PC3(xmldata)),"%"),ylab=paste0("PC4  ",get_var(find_PC4(xmldata)),"%"),
     type="p",cex=0.8,
     col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.2, "npc"), 
          y = unit(0.05, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.8, "npc"), 
          y = unit(0.05, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.05, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.05, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))

# plot PC1PC5
plot(DF_piggies$PC1,DF_piggies$PC5,
     main="guppy on all - change with time",
     xlab=paste0("PC1  ",get_var(find_PC1(xmldata)),"%"),ylab=paste0("PC5  ",get_var(find_PC5(xmldata)),"%"),
     type="p",cex=0.8,
     col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.2, "npc"), 
          y = unit(0.05, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.8, "npc"), 
          y = unit(0.05, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC5
grid.text(PC_down(find_PC5(xmldata)), x = unit(0.05, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC5(xmldata)), x = unit(0.05, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# control
par(mfrow=c(3,2), mai = c(0.4, 0.4, 0.4, 0.4))
plot(DF_piggies$PC1[DF_piggies$Cohort=="Control"],
     DF_piggies$PC2[DF_piggies$Cohort=="Control"],
     main="PC1 PC2 Control",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Control"],
     DF_piggies$PC4[DF_piggies$Cohort=="Control"],
     main="PC3 PC4 Control",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="D-scour"],
     DF_piggies$PC2[DF_piggies$Cohort=="D-scour"],
     main="PC1 PC2 D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="D-scour"],
     DF_piggies$PC4[DF_piggies$Cohort=="D-scour"],
     main="PC3 PC4 D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="ColiGuard"],
     DF_piggies$PC2[DF_piggies$Cohort=="ColiGuard"],
     main="PC1 PC2 ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="ColiGuard"],
     DF_piggies$PC4[DF_piggies$Cohort=="ColiGuard"],
     main="PC3 PC4 ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
par(mfrow=c(3,2), mai = c(0.4, 0.4, 0.4, 0.4))
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin"],
     main="PC1 PC2 Neomycin",
     xlab="PC1",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin"],
     DF_piggies$PC4[DF_piggies$Cohort=="Neomycin"],
     main="PC3 PC4 Neomycin",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin+D-scour"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin+D-scour"],
     main="PC1 PC2 Neomycin+D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin+D-scour"],
     DF_piggies$PC4[DF_piggies$Cohort=="Neomycin+D-scour"],
     main="PC3 PC4 Neomycin+D-scour",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     main="PC1 PC2 Neomycin+ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     DF_piggies$PC4[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     main="PC3 PC4 Neomycin+ColiGuard",
     xlab="",ylab="",
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
dev.off()

##############################
##############################

# piggies (guppied by time point)

DF_piggies_time
xml_data <- piggies 

# re-order cohort 
DF_piggies_time$Cohort <- factor(DF_piggies_time$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"
                       ))

PC1PC2_DF_piggies_time_plots <- 
  DF_piggies_time %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC3PC4_DF_piggies_time_plots <- 
  DF_piggies_time %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC1PC5_DF_piggies_time_plots <- 
  DF_piggies_time %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
  

pdf("piggies_time.pdf")
PC1PC2_DF_piggies_time_plots$plots[[1]]
# xml data to add 
# figuring out a way to do it in a cleaner way 
PC3PC4_DF_piggies_time_plots$plots
PC1PC5_DF_piggies_time_plots$plots
dev.off()


##############################
##############################


# piggies (guppied by time point)

groupA

# re-order cohort 
groupA$Cohort <- factor(groupA$Cohort, 
                        levels=c("Control", 
                                 "D-scour", 
                                 "ColiGuard",
                                 "Neomycin",
                                 "Neomycin+D-scour",
                                 "Neomycin+ColiGuard"
                        ))

PC1PC2_groupA_plots <- 
  groupA %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC3PC4_groupA_plots <- 
  groupA %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC1PC5_groupA_plots <- 
  groupA %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))


pdf("groupA_time.pdf")
PC1PC2_groupA_plots$plots
PC3PC4_groupA_plots$plots
PC1PC5_groupA_plots$plots
dev.off()



##############################
##############################


# piggies (guppied by time point)

groupB

# re-order cohort 
groupB$Cohort <- factor(groupB$Cohort, 
                        levels=c("Control", 
                                 "D-scour", 
                                 "ColiGuard",
                                 "Neomycin",
                                 "Neomycin+D-scour",
                                 "Neomycin+ColiGuard"
                        ))

PC1PC2_groupB_plots <- 
  groupB %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC3PC4_groupB_plots <- 
  groupB %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC1PC5_groupB_plots <- 
  groupB %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point() + 
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))


pdf("groupB_time.pdf")
PC1PC2_groupB_plots$plots
PC3PC4_groupB_plots$plots
PC1PC5_groupB_plots$plots
dev.off()




