#run in R-3.3.1
#a few lines of code from methylCtools bcall2beta by Hovestadt et al. 2014 were retained 
#CpG position handling and coverage and Beta calculations by Katarzyna Sikora
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

#methylation table
mpath<-commandArgs(trailingOnly=TRUE)[2]
mshort<-basename(mpath)
message(paste0("processing ",mpath))

## filtering thresholds
a.snp <- 0.25     # threshold for SNP filtereing (allelic frequency of illegal base)
a.cov <- 10        # min coverage on both strands (KS)

#read in methylation data; use data.table
require(data.table)

dm<-fread(mshort,header=FALSE,quote="",sep="\t")
colnames(dm)<-c("chr", "start", "end", "Beta","M", "U")
dm
nori<-nrow(dm)
message(paste0(nori," CpGs were extracted"))

dm$Cov<-rowSums(dm[,c("M","U")])
head(dm)

##filter for low coverage
message("filtering low coverage calls ..", appendLF=FALSE)
dm$Beta[dm$Cov < a.cov] <- NA
message(" done")
tot<-nori
filt.cov<-sum(is.na(dm$Beta))
pct.cov<-filt.cov/tot*100

message(sprintf("%i (%.2f%%) CpGs were filtered out to coverage less than %i (sum on both strands)",filt.cov,pct.cov,a.cov))

dm$ms<-paste(dm$chr,dm$start,sep="_")

## write output files
message("writing output files ..", appendLF=FALSE)
write.table(dm,sep="\t", row.names=FALSE, quote=FALSE, file=paste0(wdir, "/", gsub("_CpG.bedGraph",".CpG.filt.bed",mshort) ))

message("done all")

