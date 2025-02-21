#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(ChIPQC)
library(yaml)
library(stringr)

bamdir<-unlist(snakemake@params[["bams"]])
peakdir<-unlist(snakemake@params[["peaks"]])
genome<-gsub("_.+","",snakemake@params[["genome"]])
wdir <- snakemake@params[["outdir"]]
blacklist<-snakemake@params[["blacklist"]]
chipdict<-snakemake@input[["chipdict"]]

setwd(wdir)

sampleSheet<-snakemake@input[["sampleSheet"]]

#take samples,marks,replicates from chipdict
yaml<-read_yaml(chipdict,as.named.list=TRUE)
ydat<-as.data.frame(do.call(rbind,lapply(yaml,as.data.frame)))
ydat$sample<-rownames(ydat)
samples<-ydat$sample


#list of supported factors
markv<-c("H3K4me1","H3K4me2","H3K4me3","H3K27ac","H3K27me3","H3K9me3","H3K36me3","H4K16ac","RAD21","CTCF")
a<-sapply(markv,function(X)grep(X,samples,ignore.case=TRUE),simplify=TRUE)
a<-a[!lapply(a,length)<1]
b<-unlist(a)
names(b)<-sub("[0-9]$","",names(b))
markv<-names(sort(b))
if(all(is.na(markv))){
  markv<-rep("All",length(samples))
}

if(all(grepl("rep",samples))){
  #regres<-regexpr("rep[0-9]?",samples)
  repv<-str_extract(samples,"rep[0-9]+")
  repv<-as.numeric(gsub("rep","",repv))
}else{
  repv<-rep(1,length(samples))
}

#check if sample sheet is NA or a file path
#first implementation: ignore sample sheet and condition and replicates
#if sample sheet is a file path: get condition and replicate information
#the check that the sample sheet file exists is taken care of by the python wrapper
if (!is.null(sampleSheet)){
  sampleinfo<-read.table(sampleSheet,header=TRUE,sep="\t",quote="")
  condv<-sampleinfo$condition[match(samples,sampleinfo$name)]
}else{
  condv<-rep("All",length(samples))
}

sampledat<-data.frame("SampleID"=samples,"Condition"=condv,"Factor"=markv,"Replicate"=repv)

#ensure that samples,bamdir and peakdir are in the same order!
    
sampledat$bamReads<-bamdir[match(samples,sub("\\.filtered.bam","",basename(bamdir)))]
message(sprintf("Provided peak files: %s", unlist(peakdir)))
sampledat$Peaks<-peakdir[match(samples,sub("\\.filtered.+","",basename(peakdir)))]

##annotation -> check for supported genome versions
message(paste0("Provided genome: ",genome))
supported_annotations<-c("hg19","hg18","mm10","mm9","ce6","dm3")
#modify genome string
if( genome %in% supported_annotations){

    annotation<-genome


   } else {stop("No matching annotation was found.")}

blist<-ifelse(file.exists(blacklist),blacklist,NULL)
QC<-ChIPQC(sampledat,annotation=annotation,mapQCth=3,blacklist=blist)
ChIPQCreport(QC,reportFolder=wdir,facet=FALSE,colourBy="Factor")

sink("sessionInfo.txt")
sessionInfo()
sink()
