.libPaths(R.home("library"))

#system(paste0('mkdir -p ',wdir)) #for debugging
wdir<-snakemake@params[["wdir"]]
if (!dir.exists(wdir)) dir.create(wdir)
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

#set working directory
in_dirs<-snakemake@input[["indirs"]]
in_files<-unlist(lapply(in_dirs,function(X)dir(X,pattern="*.loom",full.names=TRUE)))
message(sprintf("analyzing files %s ",in_files))
samples<-snakemake@params[["samples"]]
message(sprintf("analyzing samples %s ",samples))

options(stringsAsFactors=FALSE,na.rm=TRUE,rgl.useNULL = TRUE)

library(dplyr)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

#####read-in data
l<-lapply(in_files,function(X)as.Seurat(ReadVelocity(X)))
names(l)<-samples
s<-MergeSeurat(x=l[[1]],y=unlist(l[[2:length(l)]]),add.cell.ids=names(l))
outfile<-file.path(wdir,basename(snakemake@output[["seurat"]]))
saveRDS(s,file=outfile)