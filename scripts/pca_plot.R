#!/usr/bin/Rscript
library(dplyr)
library(plotrix)

pcadat<-read.table("results/new.proj",header=T,stringsAsFactors=F)
mdat<-read.table("results/merged.tsv",sep='\t',fill=T,header=T,stringsAsFactors=F)
boggo<-inner_join(pcadat,mdat)

pdf("edge_pca.pdf")
plot(boggo$pc1,boggo$pc3,main="phylosift edge PCA on pigs",xlab="PC1",ylab="PC3",pch=NA,type="p")
cohorts=unique(sort(boggo$Cohort))
for(coho in 1:length(cohorts)){
  points(boggo$pc1[boggo$Cohort==cohorts[coho]],boggo$pc3[boggo$Cohort==cohorts[coho]],col=coho,pch=1)
}
legend("topleft",legend=cohorts, fill=palette(), cex=0.8, ncol=2)
dev.off()

#2*ylen/3
color_legend <- function(x, y, xlen, ylen, main, tiks, colors){
    text(x, y+.6, main, adj=c(0,0), cex=1.3)
    color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.8)
}

rbow <- rainbow(60, end=0.7, alpha=0.7)

pdf("pca_PC1_PC2_vs_sample_date.pdf")
plot(boggo$pc1[boggo$Cohort!="Mothers"&boggo$Cohort!="NegativeControl"],boggo$pc2[boggo$Cohort!="Mothers"&boggo$Cohort!="NegativeControl"],main="beta diversity (phylosift edge PCA)",xlab="PC1",ylab="PC2",type="p",col=rbow[as.Date(boggo$X.collection_date[boggo$Cohort!="Mothers"&boggo$Cohort!="NegativeControl"])-as.Date("2017-01-29 00:00:00")])

legvec <- c(0,15,30,45,60)
color_legend( -2.9, 4, 3.5, 1.5, "trial days:", legvec, rbow)
dev.off()
