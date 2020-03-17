######################################################################################################

# rationale: 

##### previous steps: 

# 1 # groups for guppy are made in guppy_group.R 

# 2 # guppy is run 

# 3 # .xml conversion to .txt:        <-  MUST RUN THIS BEFORE RUNNING THIS SCRIPT 

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
library(openxlsx)
library(gridExtra)
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
##############################

# run guppy_XML_process.R to get simplified df
# tried to load it with read.csv, read_csv, read.csv2, would not keep the same format!!!!! grrrrrrr

##############################
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
jplace_df_final$group <- gsub('piggies_CTRLNEO', 'groupC', jplace_df_final$group)
jplace_df_final$group <- gsub('piggies_NEONEOD', 'groupD', jplace_df_final$group)
jplace_df_final$group <- gsub('piggies_NEONEOC', 'groupE', jplace_df_final$group)
jplace_df_final$group <- gsub('piggies_CTRLDs', 'groupF', jplace_df_final$group)
jplace_df_final$group <- gsub('piggies_CTRLC', 'groupG', jplace_df_final$group)
jplace_df_final$group <- gsub('^piggies$', 'piggies_all', jplace_df_final$group)
# note to self: ^ and $ for the exact (and entire) string match 

jplace_df_final <- cSplit(jplace_df_final, "group","_")
jplace_df_final$group_1 <- gsub('pos', 'positive_controls', jplace_df_final$group_1)
jplace_df_final$group_2 <- gsub('controls', 'all_replicates', jplace_df_final$group_2)
jplace_df_final <- setnames(jplace_df_final, old = c('group_1','group_2'), new = c('sample_type','guppied_date'))

head(jplace_df_final)
unique(jplace_df_final$sample_type)
unique(jplace_df_final$guppied_date)

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
                      levels=c("positive_controls","piggies","groupA","groupB",
                               "groupC","groupD","groupE","groupF","groupG"))

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

groupC <- as.data.frame(multi_coggo$groupC)
groupC_Ja31 <- groupC %>%
  filter(guppied_date=="Ja31")
groupC_Fe7 <- groupC %>%
  filter(guppied_date=="Fe7")
groupC_Fe14 <- groupC %>%
  filter(guppied_date=="Fe14")
groupC_Fe21 <- groupC %>%
  filter(guppied_date=="Fe21")
groupC_Fe28 <- groupC %>%
  filter(guppied_date=="Fe28")
groupC_Ma3 <- groupC %>%
  filter(guppied_date=="Ma3")
groupC <- rbind(
  groupC_Ja31,
  groupC_Fe7,
  groupC_Fe14,
  groupC_Fe21,
  groupC_Fe28,
  groupC_Ma3)

##############################
##############################

groupD <- as.data.frame(multi_coggo$groupD)
groupD_Ja31 <- groupD %>%
  filter(guppied_date=="Ja31")
groupD_Fe7 <- groupD %>%
  filter(guppied_date=="Fe7")
groupD_Fe14 <- groupD %>%
  filter(guppied_date=="Fe14")
groupD_Fe21 <- groupD %>%
  filter(guppied_date=="Fe21")
groupD_Fe28 <- groupD %>%
  filter(guppied_date=="Fe28")
groupD_Ma3 <- groupD %>%
  filter(guppied_date=="Ma3")
groupD <- rbind(
  groupD_Ja31,
  groupD_Fe7,
  groupD_Fe14,
  groupD_Fe21,
  groupD_Fe28,
  groupD_Ma3)

##############################
##############################

groupE <- as.data.frame(multi_coggo$groupE)
groupE_Ja31 <- groupE %>%
  filter(guppied_date=="Ja31")
groupE_Fe7 <- groupE %>%
  filter(guppied_date=="Fe7")
groupE_Fe14 <- groupE %>%
  filter(guppied_date=="Fe14")
groupE_Fe21 <- groupE %>%
  filter(guppied_date=="Fe21")
groupE_Fe28 <- groupE %>%
  filter(guppied_date=="Fe28")
groupE_Ma3 <- groupE %>%
  filter(guppied_date=="Ma3")
groupE <- rbind(
  groupE_Ja31,
  groupE_Fe7,
  groupE_Fe14,
  groupE_Fe21,
  groupE_Fe28,
  groupE_Ma3)

##############################
##############################


groupF <- as.data.frame(multi_coggo$groupF)
groupF_Ja31 <- groupF %>%
  filter(guppied_date=="Ja31")
groupF_Fe7 <- groupF %>%
  filter(guppied_date=="Fe7")
groupF_Fe14 <- groupF %>%
  filter(guppied_date=="Fe14")
groupF_Fe21 <- groupF %>%
  filter(guppied_date=="Fe21")
groupF_Fe28 <- groupF %>%
  filter(guppied_date=="Fe28")
groupF_Ma3 <- groupF %>%
  filter(guppied_date=="Ma3")
groupF <- rbind(
  groupF_Ja31,
  groupF_Fe7,
  groupF_Fe14,
  groupF_Fe21,
  groupF_Fe28,
  groupF_Ma3)

##############################
##############################


groupG <- as.data.frame(multi_coggo$groupG)
groupG_Ja31 <- groupG %>%
  filter(guppied_date=="Ja31")
groupG_Fe7 <- groupG %>%
  filter(guppied_date=="Fe7")
groupG_Fe14 <- groupG %>%
  filter(guppied_date=="Fe14")
groupG_Fe21 <- groupG %>%
  filter(guppied_date=="Fe21")
groupG_Fe28 <- groupG %>%
  filter(guppied_date=="Fe28")
groupG_Ma3 <- groupG %>%
  filter(guppied_date=="Ma3")
groupG <- rbind(
  groupG_Ja31,
  groupG_Fe7,
  groupG_Fe14,
  groupG_Fe21,
  groupG_Fe28,
  groupG_Ma3)

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
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))
PC3PC4_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC3,y=PC4,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))
PC1PC5_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC5,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC5 (",get_var(find_PC5(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))


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

# legend 
color_legend <- function(x, y, xlen, ylen, main, tiks, colors){
  text(x, y+.55, main, adj=c(0,0), cex=1.8)
  color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.7)
}
rbow <- rainbow(40, end=0.7, alpha=0.7)
legvec <- c(0,10,20,30,40)


pdf("time_beta.pdf")
# PC1PC2
par(oma=c(6,6,6,6)) # all sides have 4 lines of space
par(mar=c(4,4,0.01,0.01))
plot(DF_piggies$PC1,DF_piggies$PC2,
     xlab=paste0("PC1  (",get_var(find_PC1(xmldata)),"%)"),
     ylab=paste0("PC2  (",get_var(find_PC2(xmldata)),"%)"),
     type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
     col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
color_legend(min(DF_piggies$PC1), max(DF_piggies$PC2)-0.6, 
             3.5, 0.9, "trial days:", legvec, rbow)
mtext(paste0(PC_down(find_PC1(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="black", outer=TRUE)  # PC1 low
mtext(paste0(PC_up(find_PC1(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="black", outer=TRUE)  # PC1 high
mtext(paste0(PC_down(find_PC2(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="black", outer=TRUE)   # PC2 low 
mtext(paste0(PC_up(find_PC2(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="black", outer=TRUE)   # PC2 high
dev.off()

# plot(DF_piggies$PC3,DF_piggies$PC4,
#      xlab=paste0("PC3  (",get_var(find_PC3(xmldata)),"%)"),
#      ylab=paste0("PC4  (",get_var(find_PC4(xmldata)),"%)"),
#      type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
#      col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
# color_legend(min(DF_piggies$PC3), max(DF_piggies$PC4)-0.25, 
#              1.8, 0.5, "trial days:", legvec, rbow)
# mtext(paste0(PC_down(find_PC3(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="blue", outer=TRUE)
# mtext(paste0(PC_up(find_PC3(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="purple", outer=TRUE)  
# mtext(paste0(PC_down(find_PC4(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="red", outer=TRUE)    
# mtext(paste0(PC_up(find_PC4(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="orange", outer=TRUE)   
# PC1PC5
# par(oma=c(6,6,6,6)) # all sides have 4 lines of space
# par(mar=c(4,4,0.01,0.01))
# plot(DF_piggies$PC1,DF_piggies$PC5,
#      xlab=paste0("PC1  (",get_var(find_PC1(xmldata)),"%)"),
#      ylab=paste0("PC5  (",get_var(find_PC5(xmldata)),"%)"),
#      type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
#      col=rbow[as.Date(DF_piggies$collection_date)-as.Date("2017-01-29 00:00:00")])
# color_legend(min(DF_piggies$PC1), max(DF_piggies$PC5)-0.15, 
#              2.7, 0.35, "trial days:", legvec, rbow)
# mtext(paste0(PC_down(find_PC1(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="blue", outer=TRUE)  
# mtext(paste0(PC_up(find_PC1(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="purple", outer=TRUE)  
# mtext(paste0(PC_down(find_PC5(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="red", outer=TRUE)    
# mtext(paste0(PC_up(find_PC5(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="orange", outer=TRUE)  

pdf("time_beta_cohorts_PC1PC2.pdf")
#################### Cohorts separately
par(oma=c(0,0,0,0)) # resetting the outer margins to default for the next plots
# control groups (Control, D-scour, ColiGuard)
par(mfrow=c(3,2), mai = c(0.4, 0.4, 0.4, 0.4))
plot(DF_piggies$PC1[DF_piggies$Cohort=="Control"],
     DF_piggies$PC2[DF_piggies$Cohort=="Control"],
     main="PC1 PC2 Control",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
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
plot(DF_piggies$PC1[DF_piggies$Cohort=="ColiGuard"],
     DF_piggies$PC2[DF_piggies$Cohort=="ColiGuard"],
     main="PC1 PC2 ColiGuard",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
# Neo groups (Neo, Neo+D, Neo+C)
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin"],
     main="PC1 PC2 Neomycin",
     xlab="PC1",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
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
plot(DF_piggies$PC1[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     DF_piggies$PC2[DF_piggies$Cohort=="Neomycin+ColiGuard"],
     main="PC1 PC2 Neomycin+ColiGuard",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies$collection_date
                               [DF_piggies$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
dev.off()

# pdf("time_beta_cohorts_PC3PC4.pdf")
# #################### Cohorts separately
# par(oma=c(0,0,0,0)) # resetting the outer margins to default for the next plots
# # control groups (Control, D-scour, ColiGuard)
# par(mfrow=c(3,2), mai = c(0.4, 0.4, 0.4, 0.4))
# plot(DF_piggies$PC3[DF_piggies$Cohort=="Control"],
#      DF_piggies$PC4[DF_piggies$Cohort=="Control"],
#      main="PC3 PC4 Control",
#      xlab="",ylab="",
#      xlim=c(-1,1.75),
#      ylim=c(-1.3,0.8),
#      cex.axis=0.8,
#      type="p",col=rbow[as.Date(DF_piggies$collection_date
#                                [DF_piggies$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
# plot(DF_piggies$PC3[DF_piggies$Cohort=="D-scour"],
#      DF_piggies$PC4[DF_piggies$Cohort=="D-scour"],
#      main="PC3 PC4 D-scour",
#      xlab="",ylab="",
#      xlim=c(-1,1.75),
#      ylim=c(-1.3,0.8),
#      cex.axis=0.8,
#      type="p",col=rbow[as.Date(DF_piggies$collection_date
#                                [DF_piggies$Cohort=="D-scour"])-as.Date("2017-01-30 00:00:00")])
# plot(DF_piggies$PC3[DF_piggies$Cohort=="ColiGuard"],
#      DF_piggies$PC4[DF_piggies$Cohort=="ColiGuard"],
#      main="PC3 PC4 ColiGuard",
#      xlab="",ylab="",
#      xlim=c(-1,1.75),
#      ylim=c(-1.3,0.8),
#      cex.axis=0.8,
#      type="p",col=rbow[as.Date(DF_piggies$collection_date
#                                [DF_piggies$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
# plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin"],
#      DF_piggies$PC4[DF_piggies$Cohort=="Neomycin"],
#      main="PC3 PC4 Neomycin",
#      xlab="",ylab="",
#      xlim=c(-1,1.75),
#      ylim=c(-1.3,0.8),
#      cex.axis=0.8,
#      type="p",col=rbow[as.Date(DF_piggies$collection_date
#                                [DF_piggies$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
# plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin+D-scour"],
#      DF_piggies$PC4[DF_piggies$Cohort=="Neomycin+D-scour"],
#      main="PC3 PC4 Neomycin+D-scour",
#      xlab="",ylab="",
#      xlim=c(-1,1.75),
#      ylim=c(-1.3,0.8),
#      cex.axis=0.8,
#      type="p",col=rbow[as.Date(DF_piggies$collection_date
#                                [DF_piggies$Cohort=="Neomycin+D-scour"])-as.Date("2017-01-30 00:00:00")])
# plot(DF_piggies$PC3[DF_piggies$Cohort=="Neomycin+ColiGuard"],
#      DF_piggies$PC4[DF_piggies$Cohort=="Neomycin+ColiGuard"],
#      main="PC3 PC4 Neomycin+ColiGuard",
#      xlab="",ylab="",
#      xlim=c(-1,1.75),
#      ylim=c(-1.3,0.8),
#      cex.axis=0.8,
#      type="p",col=rbow[as.Date(DF_piggies$collection_date
#                                [DF_piggies$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
# dev.off()



# TIME BETA DENSITIES 

mytheme <- theme(legend.position="none",
                 axis.text.x=element_text(size=4),
                 axis.title.x=element_text(size=6),
                 axis.text.y=element_text(size=4),
                 axis.title.y=element_text(size=5))

p1 <- DF_piggies %>%
  filter(collection_date=="2017-01-31"|
           collection_date=="2017-02-07"|
           collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|
           collection_date=="2017-02-28"|
           collection_date=="2017-03-03") %>%
  ggplot(., aes(x = PC1, fill = collection_date)) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies$PC1),max(DF_piggies$PC1))  +
  theme_bw()+
  mytheme +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))
g1.1 <- text_grob(paste0(PC_down(find_PC1(xmldata))),size=4,lineheight = 1)
g1.2 <- text_grob(paste0(PC_up(find_PC1(xmldata))),size=4,lineheight = 1)

p2 <- DF_piggies %>%
  filter(collection_date=="2017-01-31"|
           collection_date=="2017-02-07"|
           collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|
           collection_date=="2017-02-28"|
           collection_date=="2017-03-03") %>%
  ggplot(., aes(x = PC2, fill = collection_date)) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies$PC2),max(DF_piggies$PC2))  +
  theme_bw()+
  mytheme+
  xlab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
g2.1 <- text_grob(paste0(PC_down(find_PC2(xmldata))),size=4,lineheight = 1)
g2.2 <- text_grob(paste0(PC_up(find_PC2(xmldata))),size=4,lineheight = 1)

p3 <- DF_piggies %>%
  filter(collection_date=="2017-01-31"|
           collection_date=="2017-02-07"|
           collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|
           collection_date=="2017-02-28"|
           collection_date=="2017-03-03") %>%
  ggplot(., aes(x = PC3, fill = collection_date)) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies$PC3),max(DF_piggies$PC3))  +
  theme_bw()+
  mytheme+
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))
g3.1 <- text_grob(paste0(PC_down(find_PC3(xmldata))),size=4,lineheight = 1)
g3.2 <- text_grob(paste0(PC_up(find_PC3(xmldata))),size=4,lineheight = 1)

p4 <- DF_piggies %>%
  filter(collection_date=="2017-01-31"|
           collection_date=="2017-02-07"|
           collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|
           collection_date=="2017-02-28"|
           collection_date=="2017-03-03") %>%
  ggplot(., aes(x = PC4, fill = collection_date)) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies$PC4),max(DF_piggies$PC4))  +
  theme_bw()+
  mytheme+
  xlab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
g4.1 <- text_grob(paste0(PC_down(find_PC4(xmldata))),size=4,lineheight = 1)
g4.2 <- text_grob(paste0(PC_up(find_PC4(xmldata))),size=4,lineheight = 1)

p5 <- DF_piggies %>%
  filter(collection_date=="2017-01-31"|
           collection_date=="2017-02-07"|
           collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|
           collection_date=="2017-02-28"|
           collection_date=="2017-03-03") %>%
  ggplot(., aes(x = PC5, fill = collection_date)) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies$PC5),max(DF_piggies$PC5)) +
  theme_bw()+
  mytheme+
  xlab(paste0("PC5 (",get_var(find_PC5(xmldata)),"%)"))
g5.1 <- text_grob(paste0(PC_down(find_PC5(xmldata))),size=4,lineheight = 1)
g5.2 <- text_grob(paste0(PC_up(find_PC5(xmldata))),size=4,lineheight = 1)



for_legend_only <- DF_piggies %>%
  filter(collection_date=="2017-01-31"|
           collection_date=="2017-02-07"|
           collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|
           collection_date=="2017-02-28"|
           collection_date=="2017-03-03") %>%
  ggplot(., aes(x = PC5, fill = collection_date)) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies$PC5),max(DF_piggies$PC5)) +
  theme(legend.position="right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
leg <- get_legend(for_legend_only)

lay <- rbind(c(1,1,1,1,2,2,2,2),
             c(1,1,1,1,2,2,2,2),
             c(6,6,7,7,8,8,9,9),
             c(3,3,3,3,4,4,4,4),
             c(3,3,3,3,4,4,4,4),
             c(10,10,11,11,12,12,13,13),
             c(5,5,5,5,16,16,16,16),
             c(5,5,5,5,16,16,16,16),
             c(14,14,15,15,16,16,16,16))

pdf("time_beta_densities.pdf", width=7,height=5)
grid.arrange(p1,p2,p3,p4,p5,
             g1.1,g1.2,
             g2.1,g2.2,
             g3.1,g3.2,
             g4.1,g4.2,
             g5.1,g5.2,
             leg,
             layout_matrix = lay)
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


# groupA (guppied by time point)



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


# groupB (guppied by time point)



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


##############################
##############################
##############################


# groupC (guppied by time point)



df <- groupC # dataframe for plots to be used 
a <- "groupC" # setting for xml data extraction (only sample_type necessary) 


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


##############################
##############################
##############################


# groupD (guppied by time point)



df <- groupD # dataframe for plots to be used 
a <- "groupD" # setting for xml data extraction (only sample_type necessary) 


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


##############################
##############################
##############################


# groupE (guppied by time point)



df <- groupE # dataframe for plots to be used 
a <- "groupE" # setting for xml data extraction (only sample_type necessary) 


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



##############################
##############################
##############################


# groupF (guppied by time point)



df <- groupF # dataframe for plots to be used 
a <- "groupF" # setting for xml data extraction (only sample_type necessary) 


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

##############################
##############################
##############################

##############################
##############################
##############################


# groupG (guppied by time point)



df <- groupG # dataframe for plots to be used 
a <- "groupG" # setting for xml data extraction (only sample_type necessary) 


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


##############################
##############################
##############################

###########################################################################################
###########################################################################################
###########################################################################################

# Checking significance and plot based on significance: 


###########################################################################################

# create workbook to add stats 

wb <- createWorkbook()

###########################################################################################


###################################
###################################


# function 1. : takes dataframe, filters less than n observations, per component, per collection date, per cohort 

myfun_pvalue_df_prep1 <- function(df1_here,n) {
  # first, join the PCs to be in one columns, so you can easily apply the t.test function to one column only 
  a <- df1_here %>%
    select(PC1,PC2,PC3,PC4,PC5,Cohort,collection_date,guppied_date,sample_type) %>%
    pivot_longer(
      cols = 1:5,
      names_to = "component",
      values_to = "value",
      values_drop_na = FALSE
    )
  
  # look at distribution - filter out if less than n observations, per component, per collection date, per cohort 
  df1 <- setDT(a)[, .(Freq = .N), by = .(collection_date,Cohort,component)]
  df1 <- df1 %>% filter(!Freq<n)
  a <- merge(a,df1)
  a <- a %>% 
    select(Cohort,collection_date,guppied_date,sample_type,component,value)
  return(a)
}

###################################
###################################

# function 2. : takes dataframe, computes t.test of PC values between Cohorts (within same date,group)

myfun_pvalue_df_prep2 <- function(df2_here) {
  
  # create an empty df to build on
  df_pval <- data.frame(
    group_2 = character(),
    group_1 = character(),
    p_value = numeric(),
    which_PC = character(),
    which_colldate = character()
  )
  
  a <- df2_here
  
  # for each of the components (PC1 to PC5) ...
  listcompos <- unique(a$component)
  
  for (compo in listcompos) {
    df <- a %>% filter(component==compo)
    
    listcoldates <- unique(df$guppied_date)
    
    for (colldate in listcoldates) {
      
      z <- pairwise.t.test(df$value[df$guppied_date==colldate], df$Cohort[df$guppied_date==colldate], p.adjust = "none")$p.value
      z <- as.data.frame(z)
      z
      z$group_2 <- rownames(z)
      rownames(z) <- NULL
      z <- z %>%
        pivot_longer(
          cols = -group_2,
          names_to = "group_1",
          values_to = "p_value",
          values_drop_na = TRUE
        ) 
      
      z <- z %>% mutate(which_colldate = colldate,
                        which_PC = compo)
      df_pval <- rbind(df_pval,z) 
    }
  }
  return(df_pval)
}

###################################
###################################

# function 3. : takes p-values, removes unnecessary comparisons 

myfun_pvalues_filtering <- function(df_pval) {
  
  # Final adjustments 
  
  df_pval <- df_pval %>%
    # join the two groups for which the p-value has been computed
    mutate(comparison=paste0(group_1,"_vs_",group_2)) %>%
    select(p_value,which_colldate,which_PC,comparison,group_1,group_2)
  
  
  # filtering to keep only meaningful comparisons 
  # to be kept: 
  meaningfulcomparisons <- c("Control_vs_ColiGuard", "ColiGuard_vs_Control",
                             "Control_vs_D-scour", "D-scour_vs_Control",
                             "Control_vs_Neomycin", "Neomycin_vs_Control",
                             "Neomycin_vs_Neomycin+D-scour", "Neomycin+D-scour_vs_Neomycin",
                             "Neomycin_vs_Neomycin+ColiGuard", "Neomycin+ColiGuard_vs_Neomycin")
  
  # eliminate useless comparisons
  df_pval <- df_pval[df_pval$comparison %in% meaningfulcomparisons,]
  
  return(df_pval)
  
}

###################################
###################################

# binding all dataframes except the one where all samples (irrespective of dates) where guppied in one run 

# DF_piggies_time,groupA,groupB,groupC,groupD,groupE,groupF,groupG

all <- rbind(DF_piggies_time, 
             groupA, 
             groupB, 
             groupC, 
             groupD,
             groupE,
             groupF,
             groupG)

all$groupsplit <- paste0(all$sample_type,"_",all$guppied_date)

# splitting into multiple dataframes (by file name)
multi_DFs <- split( all , f = all$groupsplit )

# prep empty df to build on 
significant <- data.frame(
  group_1 = character(),
  group_2 = character(),
  p_value = numeric(),
  which_colldate = character(),
  which_PC = character(),
  comparison = character(),
  pval.adj = numeric(),
  groupsplit = character(),
  stringsAsFactors = FALSE
)


###################################
###################################


# runs single dataframe (one df : one guppy run) through functions; returns all p-values 


for (singl_DF in multi_DFs) {
  
  # function 1 
  a1 <- myfun_pvalue_df_prep1(singl_DF,2)
  
  # function 2
  a2 <- myfun_pvalue_df_prep2(a1)
  
  # function 3
  a3 <- myfun_pvalues_filtering(a2)
  
  # convert to class dataframe
  a4 <- as.data.frame(a3)
  
  # assign name of dataframe to dataframe (useful for later rbinding)
  a5 <- a4 %>% 
    mutate(groupsplit = singl_DF$groupsplit[1])
  
  # rbind all
  significant <- rbind(
    significant,
    a5)
  
}

###################################
###################################

# save stats
addWorksheet(wb, "p_values")
writeData(wb, sheet = "p_values", significant, rowNames = FALSE)


###################################
###################################


# Post-hoc correction: 
final <- significant %>% 
  group_by(groupsplit) %>%
  add_tally() %>%
  mutate(threshold = 0.05/n) %>%
  filter(p_value<threshold) 
final


###################################
###################################

# save stats
addWorksheet(wb, "p_adj")
writeData(wb, sheet = "p_adj", final, rowNames = FALSE)


###################################
###################################

final <- cSplit(final, "groupsplit","_")
colnames(final)[colnames(final)=="groupsplit_1"] <- "dataframe"
colnames(final)[colnames(final)=="groupsplit_2"] <- "guppied_date"


# dummy df to associate guppied_date with collection_date
dummy <- data.frame(guppied_date = as.character(c("Ja31","Fe7","Fe14","Fe21","Fe28","Ma3")),
                    collection_date = as.character(c("2017-01-31","2017-02-07","2017-02-14","2017-02-21","2017-02-28","2017-03-03")))



# adding collection_date to dataframe
df <- inner_join(final,dummy) 



# string replacement (as we haven't included DF_piggies (single guppy run), 
# for coherence we need to specify that this data comes from DF_piggies_time (multiple guppy runs))
df$dataframe <- gsub("piggies","DF_piggies_time",df$dataframe)



############

# XML data extract

# taking all along except the guppy run with all the samples from all time points 
simplified2 <- simplified %>%
  filter(!guppied_date == "all")  
# now the only piggies are from the DF_piggies_time guppy runs
simplified2$sample_type <- gsub("piggies","DF_piggies_time",simplified2$sample_type)


###################################
###################################



mytheme <- theme(legend.position = "none",
                 axis.text.x=element_text(size=3),
                 axis.title.x=element_text(size=5),
                 axis.text.y=element_text(size=3),
                 axis.title.y=element_text(size=4),
                 plot.title = element_text(size = 6, face = "bold"))

# Plots only statistically significant observations,
# reporting xml data and dataframe (guppy run) it comes from

mygrobs <- vector('list', nrow(df))
pdf("guppy_sign_cohorts.pdf", onefile = TRUE)
for (A in rownames(df)) {
  A <- as.numeric(A) # dataframes rownames must be taken as numeric
  
  # subsetting of original dataframe based on what is statistically significant (rows of df)
  pp <- eval(as.name(paste(df$dataframe[A])))  %>%
    filter(collection_date==as.character(df$collection_date[A])) %>%
    filter(Cohort==as.character(df$group_1[A])|Cohort==as.character(df$group_2[A])) %>%
    select(PC1,PC2,PC3,PC4,PC5,Cohort,collection_date,guppied_date,sample_type) %>%
    pivot_longer(
      cols = 1:5,
      names_to = "component",
      values_to = "value",
      values_drop_na = FALSE
    ) %>%
    filter(component==as.name(paste(df$which_PC[A])))
  
  pp$Cohort <- factor(pp$Cohort, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))
  
  # save some parameters to report on plot
  title1 <- unique(pp$guppied_date)
  title2 <- unique(pp$sample_type)
  PC_lab <- unique(pp$component)
  
  # add xml data 
  # subsetting the xml data based on the significant observations dataframe
  a <- df$dataframe[A]
  b <- df$guppied_date[A]
  c <- df$which_PC[A]
  xmldata <- simplified2 %>%
    filter(sample_type==a) %>%
    filter(guppied_date==b) %>%
    filter(component==c)
  
  # build plot 
  p <- ggplot(pp, aes(x=value, fill=Cohort)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'stack')+
    ggtitle(paste0(title1,"_",title2)) +
    theme_bw()+
    mytheme+
    xlab(paste0(PC_lab," (",get_var(xmldata),"%)"))+
    scale_fill_discrete(drop=FALSE)
  
  g1 <- text_grob(paste0(PC_down(xmldata)),size=3,lineheight = 1)
  g2 <- text_grob(paste0(PC_up(xmldata)),size=3,lineheight = 1)
  
  lay <- rbind(c(1,1,1,1,1),
               c(1,1,1,1,1),
               c(2,2,2,3,3))
  
  grid.arrange(p,g1,g2, layout_matrix = lay)
  
  mygrobs[[A]]  <- grid.arrange(p,g1,g2, layout_matrix = lay)
}
dev.off()

pp <- ggplot(DF_piggies, aes(x=PC1, fill=Cohort)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')+
  theme_bw()+
  mytheme+
  theme(legend.position="right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))+
  scale_fill_discrete(drop=FALSE)
leg <- get_legend(pp)

# a selection of plots from the figure generated above 
pdf("guppy_sign_cohorts_selection.pdf")
lay <- rbind(c(1,2),
             c(3,4),
             c(5,6),
             c(7,8))
grid.arrange(mygrobs[[1]]
             ,mygrobs[[4]],
             mygrobs[[5]],
             mygrobs[[6]],
             mygrobs[[7]],
             mygrobs[[8]],
             mygrobs[[9]],
             leg,
             layout_matrix = lay)
dev.off()

################################################################################################
###########################################################################################
###########################################################################################



# save stats in workbook
saveWorkbook(wb, "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/stats_guppy.xlsx", overwrite=TRUE)

