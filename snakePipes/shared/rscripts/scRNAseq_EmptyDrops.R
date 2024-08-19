.libPaths(R.home("library"))

#system(paste0('mkdir -p ',wdir)) #for debugging
wdir<-snakemake@params[["wdir"]]
if (!dir.exists(wdir)) dir.create(wdir)
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

in_dirs<-snakemake@params[["indirs"]]
message(sprintf("analyzing folders %s ",unlist(in_dirs)))

get_samples<-function(folder){
    z<-strsplit(folder, split="\\/")[[1]]
    return(z[length(z)-3])
}

samples<-unlist(lapply(in_dirs,function(X)get_samples(X)))


library(DropletUtils)
library(dplyr)
library(Seurat)
library(data.table)

filter_empty_cells<-function(folder,sample){


    sce <- read10xCounts(folder,version="3",type="sparse",col.names=TRUE)
    br.out <- barcodeRanks(counts(sce))

    png(paste0(sample,".barcoderank.png"))
    plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col="red")

    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
        legend=c("knee", "inflection"))
    dev.off()


    set.seed(100)
    e.out <- emptyDrops(counts(sce))
    e.out

    is.cell <- e.out$FDR <= 0.01
    sum(is.cell, na.rm=TRUE) 

    table(Limited=e.out$Limited, Significant=is.cell)

    is.cell2<-is.cell
    is.cell2[is.na(is.cell2)]<-FALSE
    sce.filt<-sce[,is.cell2]

    png(paste0(sample,".logProb.isCell.png"))
    plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
        xlab="Total UMI count", ylab="-Log Probability")
    dev.off()

    ##fix rownames and make seurat object
    expdat <- assay(sce.filt, "counts")
    rownames(expdat)<-gsub("-PAR-Y","",rownames(expdat))
    features<-fread(file.path(folder,"features.tsv"),header=FALSE,sep="\t",quote="")
    rn<-features$V2[match(rownames(expdat),features$V1)]
    rn2<-make.unique(rn)
    rownames(expdat)<-rn2

    seuset<-CreateSeuratObject(counts = expdat,project=sample)

    return(seuset)

}

if(length(in_dirs)>1){
    l<-mapply(SIMPLIFY=FALSE, function(X,Y) filter_empty_cells(X,Y),X=in_dirs,Y=samples)
    names(l)<-samples
    l<-l[lapply(l,length)>0]
    if(length(l)>1){
        s<-merge(x=l[[1]],y=unlist(l[2:length(l)]),add.cell.ids=names(l))
    }else{
        s<-l[[1]]
    }
}else{
    s<-filter_empty_cells(in_dirs,samples)
}

outfile<-file.path(wdir,basename(snakemake@output[["seurat"]]))
saveRDS(s,file=outfile)

sink("sessionInfo.txt")
sessionInfo()
sink()

