#!/usr/bin/env Rscript

rpkm<-read.table("results/first.rpkm")

cons<-unique(rpkm$V2[rpkm$V1=="con"])
neos<-unique(rpkm$V2[rpkm$V1=="neo"])

concount<-vector(length=length(cons))
neocount<-vector(length=length(neos))
for(con in 1:length(cons)){
  concount[con]<-sum(rpkm$V5[rpkm$V2==cons[con]])
}
for(neo in 1:length(neos)){
  neocount[neo]<-sum(rpkm$V5[rpkm$V2==neos[neo]])
}

pdf("first_rpkm.pdf",width=4,height=6)
boxplot(concount,neocount,ylim=c(0,3300),main="all aminoglycoside resistance genes,\nestimated RPKM",names=c("control","neomycin"),sub="trial day 7, immediately post-treatment",ylab="estimated RPKM")
dev.off()
