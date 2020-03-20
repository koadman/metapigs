
# CONTENT: 

# mock community 
# mock community expected vs obtained 
# protexin (D-scour)
# coliguard
# all plots
# contamination mock community
# contamination protexin (D-scour)
# contamination coliguard


library(readr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(pheatmap)

setwd("/Users/12705859/Desktop/metapigs_base/pos.controls_reads_metaphlan/")


My_Theme = theme(
  axis.text.x = element_text(size=11),
  axis.title.x = element_text(size=9), 
  axis.text.y = element_text(size = 9), 
  axis.title.y = element_text(size = 11)) 


mock <- read.csv("mock_communities_merged_abundance_table_species.txt",
                                                              "\t", header = TRUE, stringsAsFactors=FALSE)

# filter out rowsums less than 
mock <- mock %>% filter(rowSums(mock[,-1]) > 0.1)

df <- mock

sum(df$plate_1_F8_S62_R1_001_profile)
colSums(df[,-1])

colRename <- function(x) {
  setNames(x, paste0("R", seq_len(ncol(x))))
}

df <- colRename(df[,-1])

a <- cbind(mock[,1], df)
colnames(a)[1] <- "taxa_assigned"

colSums(a[,-1])

df2 <- a %>%
  pivot_longer(
    cols = starts_with("R"),
    names_to = "sample",
    #names_prefix = "s",
    values_to = "value",
    values_drop_na = TRUE
  )

# nice barplot 
mock_community_plot <- ggplot(data=df2, aes(x=sample, y=value, fill=taxa_assigned)) +
  geom_bar(stat="identity") + 
  ggtitle("Mock community") +
  ylab("read abundance (%)") +
  My_Theme+
  labs(fill = "assigned taxa")
mock_community_plot

# eliminate failed sample 
df <- df[,-8] 

#mean per row (each row is a member of the samples)
df$mean <- rowMeans(df, na.rm=TRUE)

# standard deviation
df$sd <- apply(df[,1:8],1,sd)

# re-insert IDs
df$taxa_assigned <- mock[,1]

capture.output(
  paste0("mock community mean read abundance of 7/8 replicates"),
  df,
  append=FALSE,
  file="pos_controls_numbers.txt"
)


###################

# Mock community: showing how the gram negatives are better represented than the gram positives 

taxa <- c("B. subtilis (+)", 
             "E. cloacae (-)", 
             "E. coli (-)", 
             "E. faecium (+)", 
             "S. aureus (+)", 
             "S. epidermidis (+)", 
             "P. aeruginosa (-)") 

gram <- c("P", "N", "N","P","P", "P", "N")

CFU <- as.numeric(c("64",
                    "64",
                    "2",
                    "64",
                    "64",
                    "128", 
                    "128"))

genome_size <- as.numeric(c("4.2",
                            "4.7",
                            "4.6",
                            "2.5", # 2.48 ? 
                            "2.8",
                            "2.7",
                            "6.3"))

read_count <- as.numeric(c("2.93",   # df$mean[1] # B. subtilis
                           "38",     # df$mean[6] # E. cloacae
                           "10.12",    # df$mean[7] # E. coli
                           "0.97",       # df$mean[4] # E. faecium
                           "9.90",        # df$mean[2] # S. aureus
                           "3.54",         #  df$mean[3] # S. epidermidis
                           "26.72"))        #  df$mean[9] # P. aeruginosa

df

mock_data <- data.frame(taxa, CFU, read_count, genome_size, gram)

mock_data <- mock_data %>%
  mutate(norm_CFU=CFU*genome_size)

mock_data <- mock_data %>%
  mutate(ratio_CFU=norm_CFU/(sum(norm_CFU)))

mock_data <- mock_data %>%
  mutate(ratio_read_count=read_count/(sum(read_count)))

sum(mock_data$ratio_read_count)

mock_data2 <- mock_data %>%
  select(taxa, ratio_CFU, ratio_read_count) 

mock_data2 <- mock_data2 %>%
  pivot_longer(
    cols = starts_with("ratio"),
    names_to = "count",
    values_to = "rank",
    values_drop_na = TRUE
  )

mock_data2$taxa <- factor(mock_data2$taxa, 
                             levels=c("P. aeruginosa (-)", 
                                      "E. cloacae (-)", 
                                      "E. coli (-)",
                                      "S. epidermidis (+)",
                                      "B. subtilis (+)",
                                      "E. faecium (+)",
                                      "S. aureus (+)"))

mock_data2$count <- gsub("ratio_CFU", "expected", mock_data2$count)
mock_data2$count <- gsub("ratio_read_count", "observed", mock_data2$count)

colnames(mock_data2)[colnames(mock_data2) == 'taxa'] <- 'taxa_assigned'
  
# here it's visible how the type of gram influences the read count 
q <- ggplot(data=mock_data2, aes(x=count, y=rank, fill=taxa_assigned)) +
  geom_bar(stat="identity", color="black", position=position_stack())+
  labs( y ="relative abundance (ratio)") +
  theme_minimal() +
  My_Theme+
  theme(legend.position="right",
        axis.title.x=element_blank())
q

ggsave("mock_exp_vs_obs.pdf",q, 
       width=3, height=3, units="in", scale=3)

capture.output(
  paste0("Caption mock community expected vs obtained"),
  paste0("Caption: Expected and observed relative abundance of mock community members. 
        Expected relative abundance is derived by CFU normalized by genome size. 
        Observed relative abundance is derived by reads mapping with MetaPhlAn2 normalized by genome size. 
        A higher ratio of gram negative observed than expected compared to gram positives "),
  append=TRUE,
  file="pos_controls_numbers.txt"
)

#####################################################


protexin <- read.csv("protexin_merged_abundance_table_species.txt",
                 "\t", header = TRUE, stringsAsFactors=FALSE)

# filter out rowsums less than 
protexin <- protexin %>% filter(rowSums(protexin[,-1]) > 0.1)

df <- protexin
df$ID

colRename <- function(x) {
  setNames(x, paste0("R", seq_len(ncol(x))))
}

df <- colRename(df[,-1])

df <- cbind(protexin[,1], df)
colnames(df)[1] <- "taxa_assigned"

colSums(df[,-1])

df2 <- df %>%
  pivot_longer(
    cols = starts_with("R"),
    names_to = "sample",
    #names_prefix = "s",
    values_to = "value",
    values_drop_na = TRUE
  )

# nice barplot 
Protexin_plot <- ggplot(data=df2, aes(x=sample, y=value, fill=taxa_assigned)) +
  geom_bar(stat="identity") + 
  ggtitle("D-scour") +
  ylab("read abundance (%)") +
  My_Theme+
  labs(fill = "taxa_assigned")
Protexin_plot

summary(df)

#mean per row (each row is a member of the samples)
df$mean <- rowMeans(df[,2:9], na.rm=TRUE)

# standard deviation
df$sd <- apply(df[,2:9],1,sd)

df$taxa_assigned
df$mean
df$sd


capture.output(
  paste0("D-scour mean read abundance"),
  df,
  append=TRUE,
  file="pos_controls_numbers.txt"
)

#####################################################

coliguard <- read.csv("coli_guard_merged_abundance_table_species.txt",
                     "\t", header = TRUE, stringsAsFactors=FALSE)

# filter out rowsums less than 
coliguard <- coliguard %>% filter(rowSums(coliguard[,-1]) > 0.05)

df <- coliguard

colRename <- function(x) {
  setNames(x, paste0("R", seq_len(ncol(x))))
}

df <- colRename(df[,-1])

df <- cbind(coliguard[,1], df)
colnames(df)[1] <- "taxa_assigned"

colSums(df[,-1])

df2 <- df %>%
  pivot_longer(
    cols = starts_with("R"),
    names_to = "sample",
    #names_prefix = "s",
    values_to = "value",
    values_drop_na = TRUE
  )

# heatmap 
ggplot(data = df2, aes(x=sample, y=taxa_assigned, fill=value)) + 
  geom_tile()

# nice barplot 
coliguard_plot <- ggplot(data=df2, aes(x=sample, y=value, fill=taxa_assigned)) +
  geom_bar(stat="identity") + 
  My_Theme+
  ggtitle("ColiGuard") +
  ylab("read abundance (%)") +
  labs(fill = "assigned taxa")
coliguard_plot

summary(df)

#mean per row (each row is a member of the samples)
df$mean <- rowMeans(df[,2:9], na.rm=TRUE)

# standard deviation
df$sd <- apply(df[,2:9],1,sd)

df$taxa_assigned
df$mean
df$sd

capture.output(
  paste0("ColiGuard mean read abundance"),
  df,
  append=TRUE,
  file="pos_controls_numbers.txt"
)


#####################################################

pdf("pos_controls_plots.pdf", onefile = TRUE)
mock_community_plot
Protexin_plot
coliguard_plot
dev.off()

all <- plot_grid(mock_community_plot, Protexin_plot, coliguard_plot,
          nrow=3, 
          labels=c("A", "B", "C"))

ggsave("pos_controls_plots_all.pdf",all, 
       width=3, height=3, units="in", scale=3)


# add caption
capture.output(
  paste0("Caption all plots"),
  paste0("Caption: Taxonomic profile of the positive controls obtained 
by mapping the reads against a 1M bacterial genomes database with MetaPhlAn2. 
Taxa of which reads are present in >0.1% are displayed. 
a. In-house made mock community; 
b. commercially available livestock probiotic D-scour; 
         c. DPI-developed probiotic ColiGuard"),
  append=TRUE,
  file="pos_controls_numbers.txt"
)


#####################################################

pigs <- read.csv("merged_abundance_table_species.txt",
                      "\t", header = TRUE, stringsAsFactors=FALSE)

# filter out rowsums less than 1 

colSums(pigs[,-1])

pigs$plate_3_F10_S270_R1_001_profile <- NULL
pigs$plate_6_F3_S502_R1_001_profile <- NULL

pigs <- pigs %>% filter(rowSums(pigs[,-1]) > 1)

df <- pigs

colRename <- function(x) {
  setNames(x, paste0("S", seq_len(ncol(x))))
}

df <- colRename(df[,-1])

df <- cbind(pigs[,1], df)
colnames(df)[1] <- "taxa_assigned"

colSums(df[,-1])

df

so <- as.matrix(df[,2:4])
names(so) <- df$taxa_assigned
so

pheatmap(so, display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 15,
         labels_row=as.character(df2$taxa_assigned))

df2 <- df %>%
  pivot_longer(
    cols = starts_with("S"),
    names_to = "sample",
    #names_prefix = "s",
    values_to = "value",
    values_drop_na = TRUE
  )

# barplot. too much diversity to display in any way (barplot nor heatmap)
pigs_plot <- ggplot(data=df2, aes(x=sample, y=value, fill=taxa_assigned)) +
  geom_bar(stat="identity") + 
  ggtitle("pigs") +
  ylab("read abundance (%)") +
  labs(fill = "taxa_assigned")
pigs_plot
# a possibility would be to get the genus level by using sep="_"
# and barplot 
pdf("three_pigs_plots.pdf", onefile = TRUE)
pigs_plot
dev.off()


#####################################

# contamination of mock community  
mock <- read.csv("mock_communities_merged_abundance_table_species.txt",
                 "\t", header = TRUE, stringsAsFactors=FALSE)


# removing the species of the mock comm (what's left now is contamination)
mock <- mock %>% filter(rowSums(mock[,-1]) < 0.1)

df <- mock
df$ID
df

colRename <- function(x) {
  setNames(x, paste0("R", seq_len(ncol(x))))
}

df <- colRename(df[,-1])

df <- cbind(mock[,1], df)
colnames(df)[1] <- "taxa_assigned"

colSums(df[,-1])

df

df2 <- df %>%
  pivot_longer(
    cols = starts_with("R"),
    names_to = "sample",
    #names_prefix = "s",
    values_to = "value",
    values_drop_na = TRUE
  )

# plot without numbers 
ggplot(data = df2, aes(x=sample, y=taxa_assigned, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral")

df2 <- df 
df2$taxa_assigned

so <- as.matrix(df2[,2:9])
names(so) <- df$taxa_assigned
so

# heat map with numbers 
pdf("contam_mock.pdf")
pheatmap(so, display_numbers = T,
         main="Mock community contamination",
         cluster_rows = F, cluster_cols = F, fontsize_number = 11,
         labels_row=as.character(df$taxa_assigned))
dev.off()

#####################################

# contamination of protexin? 
protexin <- read.csv("protexin_merged_abundance_table_species.txt",
                     "\t", header = TRUE, stringsAsFactors=FALSE)

df <- protexin
df$ID
df
protexin$ID

colRename <- function(x) {
  setNames(x, paste0("R", seq_len(ncol(x))))
}

df <- colRename(df[,-1])

df <- cbind(protexin[,1], df)
colnames(df)[1] <- "taxa_assigned"

# removing the species of the probiotic (what's left now is contamination)
df <- df[!grepl("Bifidobacterium_bifidum", df$taxa_assigned),]
df <- df[!grepl("Enterococcus_faecium", df$taxa_assigned),]
df <- df[!grepl("Lactobacillus_helveticus", df$taxa_assigned),]
df <- df[!grepl("Lactobacillus_delbrueckii", df$taxa_assigned),]
df <- df[!grepl("Lactobacillus_plantarum", df$taxa_assigned),]
df <- df[!grepl("Lactobacillus_rhamnosus", df$taxa_assigned),]
df <- df[!grepl("Streptococcus_thermophilus", df$taxa_assigned),]

colSums(df[,-1])

df2 <- df %>%
  pivot_longer(
    cols = starts_with("R"),
    names_to = "sample",
    #names_prefix = "s",
    values_to = "value",
    values_drop_na = TRUE
  )

# heatmap without numbers
ggplot(data = df2, aes(x=sample, y=taxa_assigned, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral")

df2 <- df 
df2$ID

so <- as.matrix(df2[,2:9])
names(so) <- df$taxa_assigned
so

# heatmap with numbers  
pdf("contam_dscour.pdf")
pheatmap(so, display_numbers = T,
         main="D-scour contamination",
         cluster_rows = F, cluster_cols = F, fontsize_number = 11,
         labels_row=as.character(df$taxa_assigned))
dev.off()

#####################################

# ColiGuard REAL contamination check: 

coliguard <- read.csv("coli_guard_merged_abundance_table_species.txt",
                      "\t", header = TRUE, stringsAsFactors=FALSE)

df <- coliguard
df$ID
df
coliguard$ID

colRename <- function(x) {
  setNames(x, paste0("R", seq_len(ncol(x))))
}

df <- colRename(df[,-1])

df <- cbind(coliguard[,1], df)
colnames(df)[1] <- "taxa_assigned"

colSums(df[,-1])


# removing the two species of the probiotic (what's left now is contamination)
df <- df[!grepl("Lactobacillus_plantarum", df$taxa_assigned),]
df <- df[!grepl("Lactobacillus_salivarius", df$taxa_assigned),]


df2 <- df %>%
  pivot_longer(
    cols = starts_with("R"),
    names_to = "sample",
    #names_prefix = "s",
    values_to = "value",
    values_drop_na = TRUE
  )

# heatmap without numbers
ggplot(data = df2, aes(x=sample, y=taxa_assigned, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral")

# heatmap 
df2 <- df

sum(df2$R5)
sum(df2$R7)

so <- as.matrix(df2[,2:9])
names(so) <- df$taxa_assigned
so


pdf("contam_coliguard.pdf")
pheatmap(so, display_numbers = T,
         main="ColiGuard contamination",
         cluster_rows = F, cluster_cols = F, fontsize_number = 11,
         labels_row=as.character(df$taxa_assigned))
dev.off()



# add caption to contamination figures
capture.output(
  paste0("Caption mock community contamination"),
  paste0("In the mock community technical replicates, 
         Lactobacillus salivarius is found in one replicate at 0.01% of the total reads"),
  paste0("Caption D-scour contamination"),
  paste0("The D-scour technical replicates contained 25 contaminants, 
of which 18 and 7 were identified at the species and at genus level, respectively. 
Contaminants were present majorly in three technical replicates (R3, R7, R8) and the most frequent contaminant (Metahobrevibacter) was present in 5 of the 8 replicates"),
  paste0("Caption ColiGuard contamination"),
  paste0("ColiGuard contained 20 contaminants, 
         of which 16 and 4 were identified at the species and at genus level, respectively. 
         Contaminants were present majorly in two technical replicates (R5, R7)."),
  append=TRUE,
  file="pos_controls_numbers.txt"
  )
