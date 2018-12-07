#!/usr/bin/Rscript
library(dplyr)

fpddat<-read.table("results/fpdalpha_div.tsv",header=T,stringsAsFactors=F)
mdat<-read.table("results/merged.tsv",sep='\t',fill=T,header=T,stringsAsFactors=F)
boggo<-inner_join(fpddat,mdat)


pdf("alpha_phylodiv.pdf",width=9,height=5)
par(mar=(c(5, 10, 4, 2) +0.1))
boxplot(boggo$phylo_entropy~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin","Neomycin+D-scour","Mothers","PosControl_ColiGuard","PosControl_Protexin","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Phylogenetic entropy",las=1)
boxplot(boggo$bwpd~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin","Neomycin+D-scour","Mothers","PosControl_ColiGuard","PosControl_Protexin","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Balance-weighted PD",las=1)
boxplot(boggo$unrooted_pd~factor(boggo$Cohort,c("Control","ColiGuard","D-scour","Neomycin+ColiGuard","Neomycin","Neomycin+D-scour","Mothers","PosControl_ColiGuard","PosControl_Protexin","MockCommunity")),horizontal=TRUE,main="Alpha diversity",xlab="Unrooted PD",las=1)
dev.off()



pdf("alpha_timeseries.pdf",width=9,height=5)
pigs=unique(boggo$isolation_source)
first=1
cohorts=unique(boggo$Cohort)
for(pig in pigs){
  if(bogsort$Cohort[bogsort$isolation_source==pig]=="Neomycin" | bogsort$Cohort[bogsort$isolation_source==pig]=="Control"){
    if(first==1){
      plot(as.Date(bogsort$X.collection_date[bogsort$isolation_source==pig]), bogsort$phylo_entropy[bogsort$isolation_source==pig],pch=NA,type="p",ylim=c(min(boggo$phylo_entropy),max(boggo$phylo_entropy)))
      first=0
    }
    lines(as.Date(bogsort$X.collection_date[bogsort$isolation_source==pig]), bogsort$phylo_entropy[bogsort$isolation_source==pig],col=which(cohorts==bogsort$Cohort[bogsort$isolation_source==pig])-4)
  }
}
dev.off()
