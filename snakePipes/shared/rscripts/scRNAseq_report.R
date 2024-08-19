.libPaths(R.home("library"))

wdir<-snakemake@params[["wdir"]]
if (!dir.exists(wdir)) dir.create(wdir)
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

rdir<-snakemake@params[["input"]]
rshort<-snakemake@params[["samples"]]

rl<-vector("list",)
for(i in seq_along(rdir)){
    tabi<-read.table(rdir[i],header=FALSE,sep=",",quote="",as.is=TRUE)
    colnames(tabi)<-c("Metric",rshort[i])
    rl[[i]]<-tabi
}

rdf <- Reduce(function(x, y, ...) merge(x, y, all = TRUE, by="Metric", sort=FALSE, ...),rl)

outf<-file.path(wdir,basename(snakemake@output[["report"]]))
write.table(rdf,outf,row.names=FALSE,quote=FALSE,sep="\t")

sink("sessionInfo.txt")
sessionInfo()
sink()
