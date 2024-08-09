.libPaths(R.home("library"))

#system(paste0('mkdir -p ',wdir)) #for debugging
wdir<-snakemake@params[["wdir"]]
if (!dir.exists(wdir)) dir.create(wdir)
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

#set working directory
in_dirs<-snakemake@params[["indirs"]]
message(sprintf("analyzing folders %s ",unlist(in_dirs)))
samples<-snakemake@params[["samples"]]
message(sprintf("analyzing samples %s ",samples))

options(stringsAsFactors=FALSE,na.rm=TRUE,rgl.useNULL = TRUE)

library(dplyr)
library(Seurat)

#####read-in data
#l<-lapply(in_dirs,function(X)Read10X(X))
#names(l)<-samples
#s<-MergeSeurat(x=l[[1]],y=unlist(l[[2:length(l)]]),add.cell.ids=names(l))
datav<-unlist(in_dirs)
names(datav)<-samples
print(datav)
expression_matrix <- Read10X(data.dir = datav)
s = CreateSeuratObject(counts = expression_matrix)
outfile<-file.path(wdir,basename(snakemake@output[["seurat"]]))
saveRDS(s,file=outfile)