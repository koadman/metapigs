#to be run before this: 
#./metapigs/scripts/add_extraction_plates_data.py 
#/Users/12705859/metapigs/source_data/pigs_samples_IDs_for_NCBI.xlsx 
#/Users/12705859/metapigs/source_data/DNA_plates.xlsx 
#/Users/12705859/metapigs/source_data/cohorts.xlsx > new_table

#input files for the script below: 
#new_table (created as above)
#rpkm_results (created with rpkm.py)

setwd("~/Desktop/rpkm_results_dir")

#remember to import the datasets this way: 
#new_table: from text (base)
#rpkm_results: from text (readr) choosing Tab as separator (if imported as (base) the number of entries for some reason will be halved )
#or import datasets using code: 
library(readr)
library(tidyr)
library(RSQLite)
new_table <- read.delim("~/Desktop/rpkm_results_dir/new_table")
rpkm_results <- read_delim("rpkm_results", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
rpkm_results0_000000000000001 <- read_delim("rpkm_results0_000000000000001", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
rpkm_results0_000000000000000000000000000001 <- read_delim("rpkm_results0_000000000000000000000000000001", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

#transform to dataframes
df1=data.frame(new_table)
df2=data.frame(rpkm_results)
df2.1=data.frame(rpkm_results0_000000000000001)
df2.2=data.frame(rpkm_results0_000000000000000000000000000001)

#give header to rpkm_results. 
names(df2)[1] <- "plate_well"
names(df2)[2] <- "gene_hit"
names(df2)[3] <- "rpkm"
names(df2.1)[1] <- "plate_well"
names(df2.1)[2] <- "gene_hit"
names(df2.1)[3] <- "rpkm"
names(df2.2)[1] <- "plate_well"
names(df2.2)[2] <- "gene_hit"
names(df2.2)[3] <- "rpkm"

#library(sqldf)

#join columns containing plate and well ID into one, without removing the original columns
df3 <- unite(df1, col = "plate_well", DNA_plate:DNA_well, sep = "_", remove = FALSE)

#merge dataframes new_table(NCBI) and rpkm_results into new dataframe, based on column called "plate_well"
df4 <- merge(df2, df3, by="plate_well")
df4.1 <- merge(df2.1, df3, by="plate_well")
df4.2 <- merge(df2.2, df3, by="plate_well")
View(df4.2)

#first replicate the column (to make sure the changing of the names in that column is working right)
df4$X.collection_date2 = df4$X.collection_date
df4.1$X.collection_date2 = df4.1$X.collection_date
df4.2$X.collection_date2 = df4.2$X.collection_date

#then change the (created) column attributes
levels(df4$X.collection_date2) <- c("01/30", "01/31", "02/01", "02/03", "02/06", "02/07", "02/08", "02/10", "02/14", "02/16", "02/17", "02/21", "02/24", "02/28", "03/03", "03/06", "03/07", "03/08", "03/09", "03/10", "08/14", "18/01/24", "NaT")
levels(df4.1$X.collection_date2) <- c("01/30", "01/31", "02/01", "02/03", "02/06", "02/07", "02/08", "02/10", "02/14", "02/16", "02/17", "02/21", "02/24", "02/28", "03/03", "03/06", "03/07", "03/08", "03/09", "03/10", "08/14", "18/01/24", "NaT")
levels(df4.2$X.collection_date2) <- c("01/30", "01/31", "02/01", "02/03", "02/06", "02/07", "02/08", "02/10", "02/14", "02/16", "02/17", "02/21", "02/24", "02/28", "03/03", "03/06", "03/07", "03/08", "03/09", "03/10", "08/14", "18/01/24", "NaT")

#normalize rpkm
df4.2$rpkm <- sqrt(df4.2$rpkm)

#df9 <- df4[, c(1, 2, 3, 22, 31, 35)]
df9.2 <- df4.2[, c(1, 2, 3, 22, 31, 35)]

#aaron code

ttester <- function(gene) {
  neo <- bobo$delta_rpkm[bobo$Cohort=="Neomycin" & bobo$gene_hit==gene]
  control <- bobo$delta_rpkm[bobo$Cohort=="Control" & bobo$gene_hit==gene]
  ttt <- 1
  if(length(neo)>2&length(control)>2){
    tttest<-t.test(neo, control)  
    ttt<-tttest$p.value
  }
  ttt
}
View(dodo_t0)

dodo_t0<-df9.2[df9.2$X.collection_date2=="01/31",c(2,3,4,5,6)]
dodo_t1<-df9.2[df9.2$X.collection_date2=="02/07",c(2,3,4,5,6)]
bobo <- merge(dodo_t0,dodo_t1, by=c("isolation_source","gene_hit","Cohort"))
#normalize rpkm.x and rpkm.y with square root
bobo$sqrt_rpkm.x <- sqrt(bobo$rpkm.x)
bobo$sqrt_rpkm.y <- sqrt(bobo$rpkm.y)
#bobo$delta_rpkm<-bobo$rpkm.y-bobo$rpkm.x
bobo$delta_rpkm<-bobo$sqrt_rpkm.y-bobo$sqrt_rpkm.x
View(bobo)


ttt<-sapply(unique(bobo$gene_hit),ttester)

summary(ttt)
View(ttt)

#pick genes for which p value <= 0.05
df10=data.frame(ttt)
names(df10)[0] <- "gene"
names(df10)[1] <- "pvalue"
df11 <- subset(df10, pvalue <= 0.05)
View(df11)


#filtering out to leave only control vs neo and pre and post neo dates
a <- df4.2 %>% filter(
  X.collection_date2 == "01/31" | X.collection_date2 == "02/07" | X.collection_date2 == "02/14",
  Cohort == "Control" | Cohort == "Neomycin" | Cohort == "D-scour" | Cohort == "ColiGuard" | Cohort == "D-scour+Neomycin" | Cohort == "ColiGuard+Neomycin",
  gene_hit == "aadA12_1_AY665771_1" | gene_hit == "aadA15_1_DQ393783_1" | gene_hit == "aadA17_1_FJ460181_1" | gene_hit == "aadA2_2_JQ364967_1" | gene_hit == "aadA3_1_AF047479_1" | gene_hit == "aminoglycosideresistanceprotein\"/protein_id=\"YP_001965793.1\"" | gene_hit == "streptomycin3''-adenylyltransferase\"/protein_id=\"CAA48308.1\""
)
View(a)

#all genes
ggplot(a, aes(x=Cohort, y=rpkm, fill=X.collection_date2)) + 
  geom_boxplot() 














#different plot per Cohort, time vs rpkm
library(reshape2)
df4.2.2.1<-melt(df4.2.2)
View(df4.2.2.1)
ggplot(df4.2.2.1,aes(x=X.collection_date2,y=value,fill=X.collection_date2))+geom_boxplot()+facet_wrap(~Cohort)


#good plot to see all cohorts; time vs rpkm 
library(ggplot2)
library(RColorBrewer)
ggplot(df4.2.2.1,aes(x=X.collection_date2,y=value,fill=X.collection_date2))+geom_boxplot(ylim=c(1,6,1))+facet_grid(~Cohort)+
  labs(x = "collection_date", y = "rpkm",color="collection_date" )+
  scale_y_continuous(limits = c(-0, 20))+
  scale_fill_brewer(palette="Purples")+
  theme_bw()+
  theme(strip.background=element_rect(fill="black"))+
  theme(strip.text=element_text(color="white", face="bold"))


#install.packages(plotly)
library(plotly)
p <- plot_ly(ggplot2::diamonds, x = ~Cohort, y = ~rpkm, color = ~X.collection_date2, type = "box") %>%
  layout(boxmode = "group")







par(mfrow = c(3,2))
#pdf("p<0.05_genes_evalue_1e-30_rpkm_vs_time.pdf")
#boxplot(df4$rpkm[df4$Cohort=="Neomycin"]~df4$X.collection_date2[df4$Cohort=="Neomycin"], main = "1e-8", las = 2)
#boxplot(df4.1$rpkm[df4.1$Cohort=="Neomycin"]~df4.1$X.collection_date2[df4.1$Cohort=="Neomycin"], main = "1e-15", las = 2)
boxplot(df4.2$rpkm[df4.2$gene_hit=="aadA12_1_AY665771_1"&df4.2$X.collection_date2=="01/31"]~df4.2$Cohort[df4.2$gene_hit=="aadA12_1_AY665771_1"&df4.2$X.collection_date2=="01/31"], main = "aadA12_1_AY665771_1", las = 2)
boxplot(df4.2$rpkm[df4.2$gene_hit=="aadA15_1_DQ393783_1"]~df4.2$X.collection_date2[df4.2$gene_hit=="aadA15_1_DQ393783_1"], main = "aadA15_1_DQ393783_1", las = 2)
boxplot(df4.2$rpkm[df4.2$gene_hit=="aadA17_1_FJ460181_1"]~df4.2$X.collection_date2[df4.2$gene_hit=="aadA17_1_FJ460181_1"], main = "aadA17_1_FJ460181_1", las = 2)
boxplot(df4.2$rpkm[df4.2$gene_hit=="aadA15_1_DQ393783_1"]~df4.2$X.collection_date2[df4.2$gene_hit=="aadA15_1_DQ393783_1"], main = "aadA15_1_DQ393783_1", las = 2)
boxplot(df4.2$rpkm[df4.2$gene_hit=="aadA15_1_DQ393783_1"]~df4.2$X.collection_date2[df4.2$gene_hit=="aadA15_1_DQ393783_1"], main = "aadA15_1_DQ393783_1", las = 2)
boxplot(df4.2$rpkm[df4.2$gene_hit=="aadA15_1_DQ393783_1"]~df4.2$X.collection_date2[df4.2$gene_hit=="aadA15_1_DQ393783_1"], main = "aadA15_1_DQ393783_1", las = 2)


boxplot(df4.2$rpkm[df4.2$gene_hit=="aadA12_1_AY665771_1"] ~df4.2$Cohort[df4.2$gene_hit=="aadA12_1_AY665771_1", subset=df4.2$X.collection_date2 == "01/31"], main = "aadA12_1_AY665771_1", las = 2)
#dev.off

#list <- as.list(unique(df4.2$gene_hit)

#list <- c("A", "B", ")

#plotting one gene rpkm vs time of Neo cohort only, with a break. 
df7 <- subset(df4, gene_hit == "aph(6)-Id_5_18676889_1")
par(mfrow = c(2,2))
cur_cnt <- 1
pdf("Neomycin_aph6.pdf")
for (i in df7){
  boxplot(df7$rpkm[df7$Cohort=="Neomycin"]~df7$X.collection_date2[df7$Cohort=="Neomycin"],
       #main = paste(i),
       xlab = "time",
       ylab = "rpkm",
       las = 2)
  cur_cnt >= 1
  if (cur_cnt > 3) {
    break
  }
}
dev.off()

#plotting one gene rpkm vs time of Neo cohort only, with a break. (0_000000000000001)
df7.1 <- subset(df4.1, gene_hit == "aph(6)-Id_5_18676889_1")
par(mfrow = c(2,2))
cur_cnt <- 1
pdf("Neomycin_aph6_0_000000000000001.pdf")
for (i in df7){
  boxplot(df7.1$rpkm[df7.1$Cohort=="Neomycin"]~df7.1$X.collection_date2[df7.1$Cohort=="Neomycin"],
          #main = paste(i),
          xlab = "time",
          ylab = "rpkm",
          las = 2)
  cur_cnt >= 1
  if (cur_cnt > 3) {
    break
  }
}
dev.off()


#plotting all genes rpkm vs time of Neo cohort only, with a break. 
par(mfrow = c(2,2))
cur_cnt <- 1
pdf("Neomycin.pdf")
for (i in df4){
  boxplot(df4$rpkm[df4$Cohort=="Neomycin"]~df4$X.collection_date2[df4$Cohort=="Neomycin"],
          #main = paste(i),
          xlab = "time",
          ylab = "rpkm",
          las = 2)
  cur_cnt >= 1
  if (cur_cnt > 3) {
    break
  }
}
dev.off()


#plotting all genes rpkm vs time of Neo cohort only, with a break (0_000000000000001)
par(mfrow = c(2,2))
cur_cnt <- 1
pdf("Neomycin_0_000000000000001.pdf")
for (i in df4.1){
  boxplot(df4.1$rpkm[df4.1$Cohort=="Neomycin"]~df4.1$X.collection_date2[df4.1$Cohort=="Neomycin"],
          #main = paste(i),
          xlab = "time",
          ylab = "rpkm",
          las = 2)
  cur_cnt >= 1
  if (cur_cnt > 3) {
    break
  }
}
dev.off()

#plotting all genes rpkm vs time of Neo cohort only, with a break (0_000000000000000000000000000001)
par(mfrow = c(2,2))
cur_cnt <- 1
pdf("Neomycin_0_000000000000000000000000000001.pdf")
for (i in df4.2){
  boxplot(df4.2$rpkm[df4.2$Cohort=="Neomycin"]~df4.2$X.collection_date2[df4.2$Cohort=="Neomycin"],
          #main = paste(i),
          xlab = "time",
          ylab = "rpkm",
          las = 2)
  cur_cnt <- cnt_cur + 1
  if (cur_cnt > 3) {
    break
  }
}
dev.off()


#it loops but it s picking everything at once, not going through each row. probably its the way the boxplot is made
par(mfrow = c(2,2))
cur_cnt <- 1
for (i in df4.2){
  plot(df4.2$X.collection_date2[df4.2$Cohort=="Neomycin"], df4.2$rpkm[df4.2$Cohort=="Neomycin"],
          #main = paste(i),
          xlab = "time",
          ylab = "rpkm",
          las = 2)
  cur_cnt <- cur_cnt + 1
  if (cur_cnt > 3) {
    break
  }
}

par(mfrow = c(3,1))
pdf("Neomycin_cohort_all_genes_rpkm_1e8_1e15_1e30.pdf")
boxplot(df4$rpkm[df4$Cohort=="Neomycin"]~df4$X.collection_date2[df4$Cohort=="Neomycin"], main = "1e-8", las = 2)
boxplot(df4.1$rpkm[df4.1$Cohort=="Neomycin"]~df4.1$X.collection_date2[df4.1$Cohort=="Neomycin"], main = "1e-15", las = 2)
boxplot(df4.2$rpkm[df4.2$Cohort=="Neomycin"]~df4.2$X.collection_date2[df4.2$Cohort=="Neomycin"], main = "1e-30", las = 2)
dev.off()

#order Cohorts on the x axis
df4$Cohort<-factor(df4$Cohort, levels=c("Mothers", "Control", "ColiGuard", "D-scour", "Neomycin", "Neomycin+ColiGuard", "Neomycin+D-scour"))
#plot rpkm by Cohort
pdf("rpkm_by_cohort.pdf")
plot(df4$Cohort, df4$rpkm, main="aminoglycoside_genes_rpkm", ylab="rpkm", las = 2)
dev.off()

#subset to Neomycin Cohort only
df6 <- subset(df4, Cohort == "Neomycin")
df7 <- subset(df6, gene_hit == "aph(6)-Id_5_18676889_1")

#subset to Neomycin only
df10 <- df9[df9$Cohort == "Neomycin", ]

#subset of useful columns (good for now, for debugging)
df9 <- df4[, c(1, 2, 3, 22, 31, 35)]

#to subset it to one specific gene:
df5 <- subset(df4, gene_hit == 'neomycinresistanceprotein"/protein_id="AAA88361.1"', select = c("rpkm", "Cohort", "X.collection_date"))

#plot one gene, setting y (rpkm) limit 
plot(df5$Cohort, df5$rpkm, main="neomycin_res_AAA88361.1_rpkm", ylab="rpkm", ylim=c(0,30))

#to look at which ones of the rpkm in df5 are really high
df5 %>% filter(rpkm > 400)

#look at how many hits depending on time 
#filter out the Neomycin cohort into new dataframe
#not working. collection dates are all NA
df7 <- subset(df4, Cohort == 'Neomycin', select = c("rpkm", "Cohort", "X.collection_date"))


#order based on date (otherwise not ordered logically)
df7$X.collection_date<-factor(df7$X.collection_date, levels=c("01/30", "01/31", "02/01", "02/03", "02/06", "02/07", "02/08", "02/10", "02/14", "02/16", "02/17", "02/21", "02/24", "02/28", "03/03", "03/06", "03/07", "03/08", "03/09", "03/10", "08/14", "18/01/24", "NaT" ))

plot(df7$X.collection_date, df7$rpkm, main="aminoglycoside_genes_rpkm_time")



