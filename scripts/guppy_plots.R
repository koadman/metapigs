######################################################################################################

# rationale: 

##### previous steps: 

# 1 # groups for guppy are made in guppy_group.R 

# 2 # guppy is run 

# 3 # .xml conversion to .txt:

# run forester.jar from command line to convert the .xml file to phyloXML - R readable format (.txt) : 
# this way: 
# java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy file.xml file.txt
# in a loop: 
# for fpath in /Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/*.xml; 
# do java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy "$fpath" 
# "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/$(basename "$fpath").txt"; 
# done

# 4 # .xml(txt) files are read in and parsed in guppy_XML_process.R

##### HERE : 

# 1 # .jplace files are read in and parsed

# 2 # bach effect removal 

# 3 # merges metadata

# 4 # principal component are plotted, where xml data complements the meaning


######################################################################################################


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
library(sva)   # this is combat , careful not to install COMBAT instead which is another thing

# manual settings 

removebatcheffect_allowed <- "yes"     # if yes, it removes the batch effect only where detected


######################################################################################################


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
    #single_DF <- unfactor(single_DF[])
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
      batch_removal=paste0("after"),
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

xmldata <- simplified %>%
  filter(sample_type=="pos") %>%
  group_split(component) 

PC1PC2_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC2,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.5,0.5,4,4),"cm"))
PC3PC4_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC3,y=PC4,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE)+
  theme(plot.margin=unit(c(0.5,1.5,1.5,1.7),"cm"))
PC1PC5_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC5,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC5 (",get_var(find_PC5(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE)+
  theme(plot.margin=unit(c(0.5,0.5,1.9,1.7),"cm"))


pdf("pos_controls.pdf")
### plot PC3PC4 
PC1PC2_pos_controls
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
### plot PC3PC4 
PC3PC4_pos_controls
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
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
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC5
grid.text(PC_down(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
dev.off()


##############################
##############################


# piggies (all time points)

a <- "piggies"

# df for plots
DF_piggies 

# df for xml data 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="all") %>%
  group_split(component) 


pdf("piggies_guppied_all.pdf")
# PC1PC2
par(oma=c(6,6,6,6)) # all sides have 4 lines of space
par(mar=c(4,4,0.01,0.01))
plot(DF_piggies$PC1,DF_piggies$PC2,
     xlab=paste0("PC1  (",get_var(find_PC1(xmldata)),"%)"),
     ylab=paste0("PC2  (",get_var(find_PC2(xmldata)),"%)"),
     type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
     col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
mtext(paste0(PC_down(find_PC1(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="blue", outer=TRUE)  # PC1 low
mtext(paste0(PC_up(find_PC1(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="purple", outer=TRUE)  # PC1 high
mtext(paste0(PC_down(find_PC2(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="red", outer=TRUE)   # PC2 low 
mtext(paste0(PC_up(find_PC2(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="orange", outer=TRUE)   # PC2 high
# PC3PC4
par(oma=c(6,6,6,6)) # all sides have 4 lines of space
par(mar=c(4,4,0.01,0.01))
plot(DF_piggies$PC3,DF_piggies$PC4,
     xlab=paste0("PC3  (",get_var(find_PC3(xmldata)),"%)"),
     ylab=paste0("PC4  (",get_var(find_PC4(xmldata)),"%)"),
     type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
     col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
mtext(paste0(PC_down(find_PC3(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="blue", outer=TRUE)
mtext(paste0(PC_up(find_PC3(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="purple", outer=TRUE)  
mtext(paste0(PC_down(find_PC4(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="red", outer=TRUE)    
mtext(paste0(PC_up(find_PC4(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="orange", outer=TRUE)   
# PC1PC5
par(oma=c(6,6,6,6)) # all sides have 4 lines of space
par(mar=c(4,4,0.01,0.01))
plot(DF_piggies$PC1,DF_piggies$PC5,
     xlab=paste0("PC1  (",get_var(find_PC1(xmldata)),"%)"),
     ylab=paste0("PC5  (",get_var(find_PC5(xmldata)),"%)"),
     type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
     col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
mtext(paste0(PC_down(find_PC1(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="blue", outer=TRUE)  
mtext(paste0(PC_up(find_PC1(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="purple", outer=TRUE)  
mtext(paste0(PC_down(find_PC5(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="red", outer=TRUE)    
mtext(paste0(PC_up(find_PC5(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="orange", outer=TRUE)  
# control groups (Control, D-scour, ColiGuard)
par(mfrow=c(3,2), mai = c(0.3, 0.3, 0.3, 0.3))
plot(DF_piggies$PC1[DF_piggies$Cohort=="Control"],
     DF_piggies$PC2[DF_piggies$Cohort=="Control"],
     main="PC1 PC2 Control",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Control"],
     DF_piggies$PC4[DF_piggies$Cohort=="Control"],
     main="PC3 PC4 Control",
     xlab="",ylab="",
     xlim=c(-1,1.75),
     ylim=c(-1.3,0.8),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="D-scour"],
     DF_piggies$PC2[DF_piggies$Cohort=="D-scour"],
     main="PC1 PC2 D-scour",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="D-scour"],
     DF_piggies$PC4[DF_piggies$Cohort=="D-scour"],
     main="PC3 PC4 D-scour",
     xlab="",ylab="",
     xlim=c(-1,1.75),
     ylim=c(-1.3,0.8),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="ColiGuard"],
     DF_piggies$PC2[DF_piggies$Cohort=="ColiGuard"],
     main="PC1 PC2 ColiGuard",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="ColiGuard"],
     DF_piggies$PC4[DF_piggies$Cohort=="ColiGuard"],
     main="PC3 PC4 ColiGuard",
     xlab="",ylab="",
     xlim=c(-1,1.75),
     ylim=c(-1.3,0.8),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
# Neo groups (Neo, Neo+D, Neo+C)
par(mfrow=c(3,2), mai = c(0.3, 0.3, 0.3, 0.3))
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin"],
     main="PC1 PC2 Neomycin",
     xlab="PC1",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin"],
     DF_piggies$PC4[DF_piggies$Cohort=="Neomycin"],
     main="PC3 PC4 Neomycin",
     xlab="",ylab="",
     xlim=c(-1,1.75),
     ylim=c(-1.3,0.8),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin+D-scour"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin+D-scour"],
     main="PC1 PC2 Neomycin+D-scour",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin+D-scour"],
     DF_piggies$PC4[DF_piggies$Cohort=="Neomycin+D-scour"],
     main="PC3 PC4 Neomycin+D-scour",
     xlab="",ylab="",
     xlim=c(-1,1.75),
     ylim=c(-1.3,0.8),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+D-scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     main="PC1 PC2 Neomycin+ColiGuard",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     DF_piggies$PC4[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     main="PC3 PC4 Neomycin+ColiGuard",
     xlab="",ylab="",
     xlim=c(-1,1.75),
     ylim=c(-1.3,0.8),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
dev.off()


##############################
##############################
##############################


# piggies (guppied by time point)



df <- DF_piggies_time # dataframe for plots to be used 
a <- "piggies" # setting for xml data extraction (only sample_type necessary) 


# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                       levels=c("Control", 
                                "D-scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-scour",
                                "Neomycin+ColiGuard"
                       ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))


pdf("piggies_guppied_by_time.pdf")
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###############################
Ja31 <- PC1PC2_plots$plots[[1]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ja31") %>%
  group_split(component) 
Ja31 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe7 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe7") %>%
  group_split(component) 
Fe7 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe14 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe14") %>%
  group_split(component) 
Fe14 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe21 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe21") %>%
  group_split(component) 
Fe21 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe28 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe28") %>%
  group_split(component) 
Fe28 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Ma3 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ma3") %>%
  group_split(component) 
Ma3 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
##############################
##############################
###############################
Ja31 <- PC3PC4_plots$plots[[1]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ja31") %>%
  group_split(component) 
Ja31 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe7 <- PC3PC4_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe7") %>%
  group_split(component) 
Fe7 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe14 <- PC3PC4_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe14") %>%
  group_split(component) 
Fe14 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe21 <- PC3PC4_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe21") %>%
  group_split(component) 
Fe21 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe28 <- PC3PC4_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe28") %>%
  group_split(component) 
Fe28 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Ma3 <- PC3PC4_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ma3") %>%
  group_split(component) 
Ma3 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
dev.off()




##############################
##############################
##############################


# piggies (guppied by time point)



df <- groupA # dataframe for plots to be used 
a <- "groupA" # setting for xml data extraction (only sample_type necessary) 


# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))


pdf(paste0(a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###############################
# Ja31 <- PC1PC2_plots$plots[[1]] 
# xmldata <- simplified %>%
#   filter(sample_type==a) %>%
#   filter(guppied_date=="Ja31") %>%
#   group_split(component) 
# Ja31 + 
#   xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
#   ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# # PC1
# grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
#           y = unit(0.1, "npc"),
#           gp = gpar(fontsize = 5))
# grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
#           y = unit(0.1, "npc"),
#           gp = gpar(fontsize = 5))
# # PC2
# grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
#           y = unit(0.25, "npc"),
#           gp = gpar(fontsize = 5))
# grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
#           y = unit(0.8, "npc"),
#           gp = gpar(fontsize = 5))
###############################
Fe7 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe7") %>%
  group_split(component) 
Fe7 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe14 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe14") %>%
  group_split(component) 
Fe14 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe21 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe21") %>%
  group_split(component) 
Fe21 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe28 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe28") %>%
  group_split(component) 
Fe28 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Ma3 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ma3") %>%
  group_split(component) 
Ma3 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
##############################
##############################
###############################
# Ja31 <- PC3PC4_plots$plots[[1]] 
# xmldata <- simplified %>%
#   filter(sample_type==a) %>%
#   filter(guppied_date=="Ja31") %>%
#   group_split(component) 
# Ja31 + 
#   xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
#   ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# # PC3
# grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
#           y = unit(0.1, "npc"),
#           gp = gpar(fontsize = 5))
# grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
#           y = unit(0.1, "npc"),
#           gp = gpar(fontsize = 5))
# # PC4
# grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
#           y = unit(0.25, "npc"),
#           gp = gpar(fontsize = 5))
# grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
#           y = unit(0.8, "npc"),
#           gp = gpar(fontsize = 5))
###############################
Fe7 <- PC3PC4_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe7") %>%
  group_split(component) 
Fe7 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe14 <- PC3PC4_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe14") %>%
  group_split(component) 
Fe14 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe21 <- PC3PC4_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe21") %>%
  group_split(component) 
Fe21 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe28 <- PC3PC4_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe28") %>%
  group_split(component) 
Fe28 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Ma3 <- PC3PC4_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ma3") %>%
  group_split(component) 
Ma3 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
dev.off()





##############################
##############################
##############################


# piggies (guppied by time point)



df <- groupB # dataframe for plots to be used 
a <- "groupB" # setting for xml data extraction (only sample_type necessary) 


# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))


pdf(paste0(a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###############################
Ja31 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ja31") %>%
  group_split(component)
Ja31 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe7 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe7") %>%
  group_split(component) 
Fe7 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe14 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe14") %>%
  group_split(component) 
Fe14 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe21 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe21") %>%
  group_split(component) 
Fe21 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe28 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe28") %>%
  group_split(component) 
Fe28 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Ma3 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ma3") %>%
  group_split(component) 
Ma3 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
##############################
##############################
###############################
Ja31 <- PC3PC4_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ja31") %>%
  group_split(component)
Ja31 +
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe7 <- PC3PC4_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe7") %>%
  group_split(component) 
Fe7 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe14 <- PC3PC4_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe14") %>%
  group_split(component) 
Fe14 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe21 <- PC3PC4_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe21") %>%
  group_split(component) 
Fe21 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Fe28 <- PC3PC4_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Fe28") %>%
  group_split(component) 
Fe28 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
Ma3 <- PC3PC4_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="Ma3") %>%
  group_split(component) 
Ma3 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################
dev.off()