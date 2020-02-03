.libPaths(R.home("library"))

#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))


options(stringsAsFactors=FALSE,na.rm=TRUE,rgl.useNULL = TRUE)

#library(Rcpp)
library(dplyr)
library(Seurat)
library(clustree)

#####read-in data
##access samples dict and STARsolo output dir
samples<-commandArgs(trailingOnly=TRUE)[2]
/data/processing/sikora/scRNA-seq.STARsolo/analysis/STARsolo/195-1/195-1.Solo.out/Gene/filtered
l<-lapply(folder_list,function(X)Read10X(X))
names(l)<-samples_dict
s<-MergeSeurat(x=l[[1]],y=unlist(l[[2:length(l)]]),add.cell.ids=names(l))
saveRDS(s,file="s.RDS")

##QC plots and filter data
