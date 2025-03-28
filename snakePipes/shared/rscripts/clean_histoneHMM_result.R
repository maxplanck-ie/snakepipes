library(GenomicRanges)

wdir <- snakemake@params[["outdir"]]
setwd(wdir)

input_peaks <- snakemake@params[["input_peaks"]]

reslist<-lapply(input_peaks,function(X)rtracklayer::import.gff(X))
names(reslist)<-gsub(".filtered.histoneHMM-regions.gff","",basename(input_peaks))
for(i in seq_along(reslist)){
  png(paste0(names(reslist)[i],"_avg_posterior.hist.png"))
  hist(as.numeric(mcols(reslist[[i]])$avg_posterior),main=names(reslist)[i],xlab="Average posterior probability")
  abline(v=0.5,col="red")
  dev.off()
}

filtlist<-lapply(reslist,function(X)X[mcols(X)$avg_posterior >= 0.5,])
for(i in seq_along(filtlist)){
  rtracklayer::export.gff3(filtlist[[i]],paste0(names(filtlist)[i],"_avgp0.5.gff"))
  a<-filtlist[[i]]
  mcols(a)$score<-as.numeric(mcols(a)$avg_posterior)
  rtracklayer::export.bed(a,paste0(names(filtlist)[i],"_avgp0.5.bed"))
}

sink("sessionInfo.txt")
sessionInfo()
sink()
