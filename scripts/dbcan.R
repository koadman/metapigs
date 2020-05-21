
library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)


setwd("~/Desktop/metapigs_dry/dbcan")
basedir = "~/Desktop/metapigs_dry/"

########################

# load dbcan_diamond data (output of mapping our-bins-AAs against CAZy db)
diamond <- read_table2("diamond.out", col_names = FALSE)

########################

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(basedir,"gtdbtk/gtdbtk_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))

########################

# load checkM assignments of the bins
checkm_all_nearly <- read_delim(paste0(basedir,"checkm/checkm_all_nearly"), 
                                "\t", escape_double = FALSE, col_types = cols(pigid = col_character()), 
                                trim_ws = TRUE)

########################

# DIAMOND data cleaning 

# get the enzyme ID clean
diamond <- cSplit(diamond,"X2","|")
diamond <- cSplit(diamond,"X2_2","_")
diamond <- cSplit(diamond,"X2_2_1"," ")

# subject/bin/contig
diamond <- cSplit(diamond,"X1","_")

diamond <- diamond %>%
  dplyr::select(X1_1,X1_2,X1_3,
                X2_1,X2_2_1_1,
                X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)

# save column with full enzyme ID 
diamond$enzyme <- diamond$X2_2_1_1

# split the enzyme ID from the enzyme number 
diamond <- diamond %>%
  separate(X2_2_1_1, 
           into = c("text", "num"), 
           sep = "(?<=[A-Za-z])(?=[0-9])")

# cols renaming
colnames(diamond) <- c("pig","bin","contig","CAZ_accession","enzymeNAME","enzymeDIGIT",
                       "perc_identity","align_length","mismatches","gap_openings",
                       "start_query","end_query","start_sample","end_sample",
                       "evalue","bitscore","enzymeID")
diamond$bin <- gsub(".fa","", diamond$bin)

head(diamond)

########################################################################

# QUALITY of diamond data : 
# how well our predicted proteins - matched against the CAZy database 

# our evalues 
summary(diamond$evalue)

# our percentages of identity 
summary(diamond$perc_identity)

# re-rder enzymes to follow the same aesthetics as Stewart et al 2019
diamond$enzymeNAME  = factor(
  diamond$enzymeNAME, levels=c("AA",
                               "CE",
                               "GH",
                               "GT",
                               "CBM",
                               "PL"))
gt_diamond$enzymeNAME <- as.character(gt_diamond$enzymeNAME)
                                                                     
# PLOT 
diamond_perc_identity_plot <- diamond %>%
  group_by(enzymeNAME) %>%
  ggplot(.,aes(enzymeNAME,perc_identity,fill=enzymeNAME))+
  ylab("Percentage identity")+
  geom_boxplot()+
  ylim(0,100)+
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Percentage identity against CAZy")

enzymes_proportion <- diamond %>% 
  group_by(enzymeNAME) %>% 
  tally() %>%
  mutate(perc = paste0(round(n/sum(n)*100,2),"%"))

RColorBrewer::display.brewer.all(n=6)


pdf("dbcan_perc_identity.pdf")
diamond_perc_identity_plot + 
  geom_text(data = enzymes_proportion,
          aes(enzymeNAME, 0, label = perc), vjust="inward",size=3)
dev.off()

########################################################################

# clean checkM all nearly data:

# filter >90 <5 bins 
checkm_all_nearly <- dplyr::filter(checkm_all_nearly, !grepl("Completeness",Completeness))
checkm_all_nearly <- subset(checkm_all_nearly, Completeness >= 90 & Contamination <= 5)
# some formatting 
checkm_all_nearly$Completeness <- as.numeric(checkm_all_nearly$Completeness)
checkm_all_nearly$Contamination <- as.numeric(checkm_all_nearly$Contamination)
# rename cols to matching colnames between dataframes to merge  
colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="pigid"] <- "pig"
colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="Bin Id"] <- "bin"
colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="Taxonomy (contained)"] <- "taxa"

checkm_clean <- checkm_all_nearly %>%
  select(pig,bin,taxa)
checkm_clean <- cSplit(checkm_clean, "taxa", sep=";")
colnames(checkm_clean) <- c("pig","bin","kingdom","phylum","class","order","family","genus","species")
head(checkm_clean)

checkm_clean$phylum <- as.character(checkm_clean$phylum)
checkm_clean$phylum[is.na(checkm_clean$phylum)] <- "Unknown"

checkm_clean$phylum <- gsub("p__","", checkm_clean$phylum)

########################################################################

# merge diamond data with checkM- nearly complete bins - taxonomic assignments 


NROW(checkm_clean)
NROW(diamond)

cm_diamond <- inner_join(checkm_clean,diamond)
NROW(cm_diamond)


cm_diamond_taxa <- cm_diamond %>%
  group_by(phylum,enzymeNAME) %>%
  tally() %>%
  mutate(Proportion = n) 

cm_diamond_prop <- cm_diamond %>%
  group_by(phylum) %>%
  tally() %>%
  mutate(proportion=n/sum(n)*100)

a <- ggplot(cm_diamond_prop, aes(y=proportion, x=phylum)) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
  theme(axis.text.x=element_text(angle=90))+
  theme_pubr()+
  ylab("Proportion od CAZymes (%)")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
b <- ggplot(cm_diamond_taxa, aes(fill=enzymeNAME, y=Proportion, x=phylum)) + 
  geom_bar(position="fill", stat="identity",colour="black")+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  theme_pubr()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x = element_blank(),
        legend.position="right",
        legend.title = element_blank())
  

# RColorBrewer::display.brewer.all(n=6)


pdf("dbcan_CAZ_CMphyla.pdf")
plot_grid(
  plot_grid(
    a + theme(legend.position = "none"),
    b + theme(legend.position = "none"),
    ncol = 1,
    align = "v"),
  plot_grid(
    get_legend(a),
    get_legend(b),
    ncol =1),
  rel_widths = c(8,2)
)
dev.off()


########################################################################

# merge diamond data with gtdbtk- bins taxonomic assignments 


NROW(gtdbtk_bins)
NROW(diamond)

gt_diamond <- inner_join(gtdbtk_bins,diamond)
NROW(gt_diamond)

gt_diamond$phylum[is.na(gt_diamond$phylum)] <- "Unknown"

gt_diamond_taxa <- gt_diamond %>%
  group_by(phylum,enzymeNAME) %>%
  tally() %>%
  mutate(Proportion = n) 

gt_diamond_prop <- gt_diamond %>%
  group_by(phylum) %>%
  tally() %>%
  mutate(proportion=n/sum(n)*100)

a_gt <- ggplot(gt_diamond_prop, aes(y=proportion, x=phylum)) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
  theme(axis.text.x=element_text(angle=90))+
  theme_pubr()+
  ylab("Proportion od CAZymes (%)")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
b_gt <- ggplot(gt_diamond_taxa, aes(fill=enzymeNAME, y=Proportion, x=phylum)) + 
  geom_bar(position="fill", stat="identity",colour="black")+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  theme_pubr()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x = element_blank(),
        legend.position="right",
        legend.title = element_blank())


pdf("dbcan_CAZ_GTphyla.pdf")
plot_grid(
  plot_grid(
    a_gt + theme(legend.position = "none"),
    b_gt + theme(legend.position = "none"),
    ncol = 1,
    align = "v"),
  plot_grid(
    get_legend(a_gt),
    get_legend(b_gt),
    ncol =1,
    rel_heights = c(3,7)),
  rel_widths = c(8,2)
)
dev.off()


########################################################################

gt_diamond <- as.data.frame(gt_diamond)


df <- gt_diamond %>%
  dplyr::select(species,bin,pig,date,phylum,enzymeNAME,enzymeID)









# normalize by lib size: 
no_reps_all_norm <- no_reps_all %>%
  group_by(cohort,pig,date) %>%
  mutate(norm_value=value/sum(value))


diamond_count <- diamond %>%
  group_by(pig,bin,enzymeID,enzymeNAME)%>%
  dplyr::summarise(enz_count=n())


df <- left_join(no_reps_all_norm,diamond_count) 


df1 <- df %>%
  mutate(enz_countBYmap=norm_value*enz_count)

df1 <- df1[!is.na(df1$enzymeID),]

df2 <- df1 %>%
  group_by(cohort,pig,date,enzymeID,enzymeNAME) %>%
  dplyr::summarise(norm_count=sum(enz_countBYmap))


# reorder dates 
df2$date  = factor(df2$date, levels=c("t0",
                                                  "t1", 
                                                  "t2",
                                                  "t3",
                                                  "t4",
                                                  "t5",
                                                  "t6",
                                                  "t7",
                                                  "t8",
                                                  "t9",
                                                  "t10",
                                      "tM"))



df2 <- as.data.frame(df2)

df_piggies <- df2 %>%
  filter(!cohort=="Mothers")

# split df by species
multiple_DFs <- split( df_piggies , f = df_piggies$enzymeID ,drop = TRUE)

# total number of unique CAZ IDs
NROW(multiple_DFs)

# empty df
all_pvalues <- data.frame(group1=character(),
                          group2=character(),
                          p_value=double(),
                          enzID=character(),
                          test=character(),
                          p_adj_method=character())


for (single_DF in multiple_DFs) {
  
  single <- as.data.frame(single_DF)
  
  
  if (NROW(unique(single$date))>10) {
    
    enzID <- single$enzymeID[1]
    
    out <- pairwise.wilcox.test(single$norm_count, single$date,
                                p.adjust.method = "bonferroni")
    
    out2 <- as.data.frame(out$p.value)
    out2$group1 <- rownames(out2)
    out2 <- out2 %>% pivot_longer(cols = -group1,names_to = "group2", values_to="p_value")
    out2 <- as.data.frame(out2[!is.na(out2$p_value),])
    out2$enzID = paste0(enzID)
    out2$test <- "pairwise.wilcox.test"
    out2$p_adj_method = "Bonferroni"
    
    all_pvalues = rbind(all_pvalues,out2)
    
  }
  
  else (print("not enough observations"))
}
 

NROW(all_pvalues)

significant <- all_pvalues %>%
  filter(p_value<0.04)

mylist <- unique(significant$enzID)


# now I will use these IDs to immediately plot interesting stuff


df_part1 <- subset(df2, (enzymeID %in% mylist)) %>% filter(!enzymeNAME=="GH")
multi_DFs_part1 <- split( df_part1 , f = df_part1$enzymeNAME ,drop = TRUE)

df_part2 <- subset(df2, (enzymeID %in% mylist)) %>% filter(enzymeNAME=="GH") %>% arrange(enzymeID)


pdf("xxxx.pdf")
for (single_DF in multi_DFs_part1) {
  
  DF <- as.data.frame(single_DF)

  print(DF %>%
    group_by(cohort,pig,date,enzymeID) %>%
    dplyr::summarize(tot=sum(norm_count)) %>%
    ggplot(., aes(date,log(tot)))+
    geom_boxplot() +
    facet_wrap(~enzymeID))
}
# need to find a smarter method to split this one in three
print(df_part2[1:25000,] %>%
        group_by(cohort,pig,date,enzymeID) %>%
        dplyr::summarize(tot=sum(norm_count)) %>%
        ggplot(., aes(date,log(tot)))+
        geom_boxplot() +
        facet_wrap(~enzymeID))
print(df_part2[25000:50000,] %>%
        group_by(cohort,pig,date,enzymeID) %>%
        dplyr::summarize(tot=sum(norm_count)) %>%
        ggplot(., aes(date,log(tot)))+
        geom_boxplot() +
        facet_wrap(~enzymeID))
print(df_part2[50000:76798,] %>%
        group_by(cohort,pig,date,enzymeID) %>%
        dplyr::summarize(tot=sum(norm_count)) %>%
        ggplot(., aes(date,log(tot)))+
        geom_boxplot() +
        facet_wrap(~enzymeID))
dev.off()




AA <- df_part1 %>%
  filter(enzymeNAME=="AA") %>%
  group_by(cohort,pig,date,enzymeID) %>%
  dplyr::summarize(tot=sum(norm_count)) 

AA %>%
  ggplot(., aes(date,log(tot)))+
  geom_boxplot() +
  facet_wrap(~enzymeID)


whos <- left_join(AA,gtdbtk_bins) %>%
  filter(enzymeID=="AA1")
unique(whos$genus)
unique(whos$species)

unique(diamond$pig)

diamond %>%
  filter(pig=="Protexin") %>%
  group_by(enzymeNAME) %>%
  tally() %>%
  mutate(prop=n/sum(n)*100)

diamond %>%
  filter(!pig=="Protexin"|
           pig=="MockCommunity"|
           pig=="ColiGuard") %>%
  group_by(enzymeNAME) %>%
  tally() %>%
  mutate(prop=n/sum(n)*100)


