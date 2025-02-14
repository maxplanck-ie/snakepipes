#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(ChIPQC)

bamdir<-snakemake@input[["bams"]]
peakdir<-snakemake@input[["peaks"]]
genome<-snakemake@params[["genome"]]
wdir <- snakemake@params[["outdir"]]
blacklist<-snakemake@params[["blacklist"]]

setwd(wdir)

sampleSheet<-snakemake@input[["sampleSheet"]]
samples<-snakemake@params[["samples"]]

#list of supported factors
markv<-c("H3K4me1","H3K4me2","H3K4me3","H3K27ac","H3K27me3","H3K9me3","H3K36me3","H4K16ac")
a<-sapply(markv,function(X)grep(X,samples,ignore.case=TRUE),simplify=TRUE)
a<-a[!lapply(a,length)<1]
b<-unlist(a)
names(b)<-sub("[0-9]$","",names(b))
markv<-names(b)

#check if sample sheet is NA or a file path
#first implementation: ignore sample sheet and condition and replicates
#if sample sheet is a file path: get condition and replicate information
#the check that the sample sheet file exists is taken care of by the python wrapper
if (sampleSheet == "" || is.null(sampleSheet)){
    sampledat<-data.frame("SampleID"=samples,"Condition"=rep("All",length(samples)),"Factor"=markv,"Replicate"=rep("All",length(samples))}

#skipped: check if samples,bamdir and peakdir are in the same ordera!
    
sampledat$bamReads<-bamdir[match(samples,sub(".filtered.bam","",basename(bamdir)))]
sampledat$Peaks<-peakdir[match(samples,gsub(".filtered.+","",basename(peakdir)))]

##annotation -> check for supported genome versions
supported_annotations<-c("hg19"="GRCh38","hg18"="GRCh37","mm10"="GRCm38","mm9"="GRCm37","ce6"="ce6","dm3"="dm3")
#modify genome string
if( any(supported_annotations[genome])){

    annotation<-supported_annotations[genome]


   } else {stop("No matching annotation was found.")}

blist<-ifelse(file.exists(blacklist),blacklist,NULL)
QC<-ChIPQC(sampledat,annotation=annotation,mapQCth=3,blacklist=blist)
ChIPQCreport(QC,reportFolder=wdir,facet=FALSE,colourBy="Factor")

sink("sessionInfo.txt")
sessionInfo()
sink()
