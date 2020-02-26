######################################################################################################

# rationale: 

# as a proof of concept: 
# run beta diversity analysis on pos controls only (all replicates) 
# and see how they cluster (do they separate? they absolutely should!)

# also check what effect the batch effect removal has on the clustering
# does it become better or worse? 


# 1 # prepares metadata files to feed guppy_epca.sh:

# 2 # run guppy_epca.sh on HPC --> output

# 3 # outputs: HPC -> local (/Users/12705859/Desktop/metapigs_base/phylosift/input_files/)
# then to R 

# 4 # remove batch effect (optional, plot with and without batch effect removal)

# 5 # merge with metadata 

# 6 # plot first 5 principal components, coloring by positive control 
# if separation is seen, the xlm file can be interrogated to determine:
# 6.1. the percentage variation explained by the principal component
# 6.2. the tree edges read abundance (+ and -) exaplined by the principal component

# let's the sanity check start! 


######################################################################################################


# 0 # set working directory & load libs

setwd("/Users/12705859/Desktop/metapigs_base/phylosift/guppy")

library(vcd)
library(summarytools)


# 1 # prepares metadata files to feed guppy_epca.sh:
# one group = one time point : 1 breed : max 2 days diff in bdays

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

# filter out pos controls, neg controls and mother samples 
mdat_sel <- mdat %>% 
  filter(Cohort=="MockCommunity"|Cohort=="PosControl_D-scour"|Cohort=="PosControl_ColiGuard") %>% 
  dplyr::select(isolation_source,DNA_plate,DNA_well)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

mdat_sel$isolation_source
# we have 9 Mock , 8 Protexin, 8 ColiGuard = 25 in total 

mdat_sel <- as.character(mdat_sel$ID)

writeLines(unlist(mdat_sel), "pos_controls_sel.txt", sep = " ")
# contains 25 file IDs in total 

###########################################################################################

# run phylosift guppy_epca.sh : location on HPC: /shared/homes/12705859/phylosift_metapigs_20200225
# this way: 
# cd /shared/homes/12705859/phylosift_metapigs_20200225
# /shared/homes/12705859/phylosift_v1.0.1/bin/guppy  epca --prefix pigpca_pos_controls_sel `cat pos_controls_sel.txt`
# output: 4 pigpca files:
# pigpca_pos_controls_sel.edgediff
# pigpca_pos_controls_sel.proj
# pigpca_pos_controls_sel.trans
# pigpca_pos_controls_sel.xml
# 
# then moved to /Users/12705859/Desktop/metapigs_base/phylosift/guppy to visualise with R 

# moved all pigpca* files (48 files) to /Users/12705859/Desktop/metapigs_base/phylosift/guppy

###########################################################################################


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
             legend.position="top",
             plot.margin=unit(c(0.3,0.3,0.3,0.3),"line"))

###########################################################################################


# upload pca output files

pcadat<-read_csv("pigpca_pos_controls_sel.proj",col_names = FALSE)


pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("PC1","PC2","PC3","PC4","PC5")
pcadat <- pcadat %>%
  dplyr::select(DNA_plate,DNA_well,PC1,PC2,PC3,PC4,PC5)


###########################################################################################

# check significance of batch effect BEFORE batch effect removal 
# checked if any batch effect still present and it seems fine
aov.out_check = aov(PC1 ~ DNA_plate, data=pcadat)   # checked for other components too
res <- TukeyHSD(aov.out_check)
aov.out_check <- as.data.frame(res$DNA_plate)
aov.out_check$`p adj`

# no need to remove batch effect as padj is approx. 1 

# pcadat is now ready to use 

###########################################################################################


# 5 # merge with metadata and average out duplicate samples -> boggo 

# merge metadata with beta div 

coggo <- merge(pcadat,mdat)
coggo <- coggo %>%
  dplyr::select(isolation_source,DNA_plate,DNA_well,PC1,PC2,PC3,PC4,PC5,sample_name)



# need to re-name Protexin to D-scour for consistency
coggo$isolation_source <- gsub("Protexin","D-scour",coggo$isolation_source)

# re-order isolation_source
coggo$isolation_source <- factor(coggo$isolation_source, 
                                 levels=c("MockCommunity", 
                                          "D-scour", 
                                          "ColiGuard"))

# 6 # plot first 4 principal components 

PC1PC2 <- coggo %>%
  ggplot(., aes(x=PC1,y=PC2,color=isolation_source))+
  geom_point(size=0.5)+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_manual(labels = c("MockCommunity", 
                                "D-scour",
                                "ColiGuard"), 
                     values = c("#aaaaaa", 
                                "#B79F00",
                                "#00BA38")) 

PC3PC4 <- coggo %>%
  ggplot(., aes(x=PC3,y=PC4,color=isolation_source))+
  geom_point(size=0.5)+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_manual(labels = c("MockCommunity", 
                                "D-scour",
                                "ColiGuard"), 
                     values = c("#aaaaaa", 
                                "#B79F00",
                                "#00BA38")) 

PC1PC5 <- coggo %>%
  ggplot(., aes(x=PC1,y=PC5,color=isolation_source))+
  geom_point(size=0.5)+
  geom_text(aes(label=isolation_source),hjust=0, vjust=0)+
  theme+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete()

all <- ggarrange(PC1PC2,PC3PC4,PC1PC5,common.legend = TRUE,
                 widths = c(1,1,1),heights=c(1,1,1))

pdf("pos_controls.pdf")
all
PC1PC2
PC3PC4
PC1PC5
dev.off()


# open pigpca_pos_controls.xml 

# gather % for each PC and 
# tree edges incr and decr with each PC

# PC1 83.6%
# PC2 14.1%
# PC3 1.5%
# PC4 0.4%
# PC5 0.2%

# 

require(XML)
data <- xmlParse("/Users/12705859/Desktop/metapigs_base/phylosift/guppy/pigpca_pos_controls_sel.xml")

xml_data <- xmlToList(data)

maybe_PC1 <- xml_data[1]

grep("color",maybe_PC1)
z <- as.data.frame(maybe_PC1)

rownames(z)
head(z)
NROW(z)
tr <- t(z)
rownames(tr)
head(tr,n=100)
NROW(tr)
View(tr)

location <- as.list(xml_data[["clade"]][["branch_length"]])
View(location)
start_time <- unlist(xml_data[["clade"]][["branch_length"]][
  names(xml_data[["clade"]][["branch_length"]]) == "start-valid-time"])



library(readr)


myfile <- basename("~/Desktop/test.txt")

my.basedir <- "~/Desktop/metapigs_base/phylosift/guppy/Archeopteryx_trees"

my.files = list.files(my.basedir,pattern=".txt")
for (textfile in my.files) {
  
  # read in file 
  my.df <- read_csv(file.path(my.basedir,textfile), col_names = FALSE)
  
  # get info from file name to add later to dataframe 
  PC_variation_explained <- round(as.numeric(gsub("<[^>]+>", "",my.df[4,])),4)*100
  myfile <- basename(textfile)
  
  # start the grepping of useful tree info 
  mylist <- grep("blue", my.df$X1)
  myend <- lapply(mylist, print)
  
  df <- data.frame(matrix(unlist(myend), nrow=length(myend), byrow=T))
  colnames(df) <- "G"
  mysel <- df %>%
    mutate(FF=G-1) %>%
    mutate(E=FF-1) %>%
    mutate(D=E-1) %>%
    mutate(C=D-1) %>%
    mutate(B=C-1) %>%
    mutate(A=B-1) %>%
    mutate(placeholder="placeholder")%>%
    pivot_longer(cols=A:G)
  
  mysel <- as.data.frame(mysel)
  row.names(mysel) <- mysel$value
  
  almostthere <- test[match(rownames(mysel), rownames(test), nomatch=0),]
  almostthere <- cSplit(almostthere, "X1","</")
  ofinterest <- almostthere$X1_1
  ofinterest <- as.data.frame(ofinterest)
  
  myprecious <- ofinterest %>% 
    filter(ofinterest != '') %>%  
    mutate(key = rep(c('taxa', 'branch_length', 'branch_width','color','red','green','blue'), n() / 7), 
           id = cumsum(key == 'taxa')) %>% 
    spread(key, ofinterest) 
  
  myprecious <- as.data.frame(lapply(myprecious, function(y) gsub("<[^>]+>", "", y)))
  
  myprecious <- myprecious %>%
    arrange(., desc(red))
  
  myprecious$PC_variation_explained <- PC_variation_explained
  myprecious$file <- myfile
  
  last.df <- rbind(
    last.df, 
    myprecious
    )
}

last.df


