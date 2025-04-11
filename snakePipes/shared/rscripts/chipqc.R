#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(ChIPQC)
library(yaml)
library(stringr)
library(purrr)


#options(MulticoreParam=MulticoreParam(workers=8))
register(MulticoreParam(8))
registered()$MulticoreParam

bamdir<-unlist(snakemake@params[["bams"]])
peakdir<-unlist(snakemake@params[["peaks"]])
genome<-gsub("_.+","",snakemake@params[["genome"]])
wdir <- snakemake@params[["outdir"]]
blacklist<-snakemake@params[["blacklist"]]
chipdict<-snakemake@input[["chipdict"]]

setwd(wdir)

spikein<-toupper(snakemake@params[["useSpikeinForNorm"]])
message(paste0("useSpikeinForNorm is set to: ",spikein))
if(spikein){
    ms<-"host"}else{ms<-"filtered"}


sampleSheet<-snakemake@input[["sampleSheet"]]

#take samples,marks,replicates from the union of narrow samples and broad samples

#yaml<-read_yaml(chipdict,as.named.list=TRUE) #not used due to buggy conversion of yaml -> list with NULL entries -> data.frame: all-NULL entries are dropped entirely
narrow_samples<-unlist(snakemake@params[["narrow_samples"]])
broad_samples<-unlist(snakemake@params[["broad_samples"]])
samples<-c(narrow_samples,broad_samples)
ydat<-data.frame("sample"=samples,"broad"=c(rep(FALSE,length(narrow_samples)),rep(TRUE,length(broad_samples))))
rownames(ydat)<-ydat$sample

ydat

#list of supported factors
markv<-c("H3K4me1","H3K4me2","H3K4me3","H3K27ac","H3K27me3","H3K9me3","H3K36me3","H4K16ac","RAD21","CTCF","MSL2","BMAL1","CLOCK")
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
    
sampledat$bamReads<-bamdir[match(samples,sub(paste0("\\.",ms,".bam"),"",basename(bamdir)))]
message(sprintf("Provided peak files: %s", unlist(peakdir)))
##for MACS2, modify input peak files: .xls -> .narrowPeak, .broadPeak
if(all(grepl("histoneHMM",peakdir))){
sampledat$Peaks<-peakdir[match(samples,sub("_avgp0.5.bed","",basename(peakdir)))]
}else{sampledat$Peaks<-peakdir[match(samples,sub("\\.filtered.+","",basename(peakdir)))]}

sampledat$PeakCaller<-"bed"
sampledat$PeakFormat<-"bed"	
if(all(grepl("MACS2",sampledat$Peaks))){
        
        #samples should be in the same order
        sampledat$Peaks[ydat$broad==TRUE]<-gsub(paste0(".",ms,".BAM_peaks.xls"),paste0(".",ms,".BAM_peaks.broadPeak"),sampledat$Peaks[ydat$broad==TRUE])
        sampledat$Peaks[ydat$broad==FALSE]<-gsub(paste0(".",ms,".BAM_peaks.xls"),paste0(".",ms,".BAM_peaks.narrowPeak"),sampledat$Peaks[ydat$broad==FALSE])
        sampledat$PeakFormat[ydat$broad==FALSE]<-"narrow"
        sampledat$PeakCaller[ydat$broad==FALSE]<-"narrow"
}

sampledat

##annotation -> check for supported genome versions
message(paste0("Provided genome: ",genome))
supported_annotations<-c("hg19","hg18","mm10","mm9","ce6","dm3")
extended_annotations<-c("GRCh38","GRCh37","GRCm38","GRCm37","ce6","dm3")
#modify genome string
if( genome %in% supported_annotations){

    annotation<-genome
} else if (genome %in% extended_annotations){
 
    annotation<-supported_annotations[grep(genome,extended_annotations)]
    
}else {stop("No matching annotation was found.")}

if(file.exists(blacklist)){
    blist<-blacklist}else{blist<-NULL}

message(paste0("Using blacklist: ",blist))
QC<-ChIPQC(sampledat,annotation=annotation,mapQCth=3,blacklist=blist)
ChIPQCreport(QC,reportFolder=".",facet=FALSE,colourBy="Factor")

sink("sessionInfo.txt")
sessionInfo()
sink()
