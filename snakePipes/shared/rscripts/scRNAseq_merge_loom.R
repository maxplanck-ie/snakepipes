sink(snakemake@log[["out"]])
.libPaths(R.home("library"))

#set working directory
in_files<-snakemake@input[["ifiles"]]
message(sprintf("analyzing files %s ",unlist(in_files)))
samples<-snakemake@params[["samples"]]
message(sprintf("analyzing samples %s ",samples))

#system(paste0('mkdir -p ',wdir)) #for debugging
wdir<-snakemake@params[["wdir"]]
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

options(stringsAsFactors=FALSE,na.rm=TRUE,rgl.useNULL = TRUE)

library(dplyr)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

#####read-in data
l<-lapply(in_files,function(X)as.Seurat(ReadVelocity(X)))
names(l)<-samples
s<-MergeSeurat(x=l[[1]],y=unlist(l[[2:length(l)]]),add.cell.ids=names(l))
outfile<-filepath(wdir,basename(snakemake@output[["seurat"]]))
saveRDS(s,file=outfile)

message('done all')
sink()
