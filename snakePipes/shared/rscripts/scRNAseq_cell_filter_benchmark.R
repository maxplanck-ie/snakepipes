#run in R3.4.0
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))


options(stringsAsFactors=FALSE,na.rm=TRUE)

require(data.table)
require(ggplot2)
require(stringr)
require(monocle)

#####read-in de-multiplexed count table
csvpath<-commandArgs(trailingOnly=TRUE)[2]

csv.dat<-fread(csvpath,header=TRUE,sep="\t")
#impute 0 where NA
csv.dat[is.na(csv.dat)]<-0

x<-csv.dat
rownames(x) <- x$GENEID
setkey(x,GENEID)
expdat <- as.data.frame(x[grep("ERCC",rownames(x),invert=TRUE),-1,with=FALSE])
rownames(expdat)<-rownames(x)[grep("ERCC",rownames(x),invert=TRUE)]
rm(x)

save(csv.dat,file="csv.dat.RData")


z<-colSums(expdat)
summary(z)

png("Expdata.ColumnSums.png")
    plot(z[order(z,decreasing=TRUE)])
    for(i in seq_along(quantile(z))){
        abline(h=quantile(z)[i],col="red")
    }
dev.off()

##make sample info

sampleInfo<-data.frame(colnames(expdat),stringsAsFactors=FALSE)
colnames(sampleInfo)<-"SampleID"
rownames(sampleInfo)<-sampleInfo$SampleID
save(sampleInfo,file="sampleInfo.RData")




#######################
## presets:
minE=2
minN=4
minT<-c(1000,1500,2000,2500,3000,3500,4000,4500,5000)

metrics.tab<-as.data.frame(minT,stringsAsFactors=FALSE)
colnames(metrics.tab)<-"minT"
metrics.tab$medGPC<-NA
metrics.tab$num_cells<-NA
metrics.tab$gene_universe<-NA

for(i in seq_along(minT)){
    minTi=minT[i]
    
####monocle
    pd <- new("AnnotatedDataFrame", data = sampleInfo)
    mono.set <- newCellDataSet(as.matrix(expdat),phenoData = pd ,lowerDetectionLimit=1,expressionFamily=negbinomial.size())
    save(mono.set,file=paste0("minT",minTi,".mono.set.RData"))
    
    pData(mono.set)$TPC <- Matrix::colSums(exprs(mono.set))
    mono.set<-mono.set[,pData(mono.set)$TPC >= minTi]
    metrics.tab$num_cells[i]<-ncol(mono.set)

    head(pData(mono.set))

    mono.set <- estimateSizeFactors(mono.set)
    mono.set <- estimateDispersions(mono.set)
    mono.set <- detectGenes(mono.set, min_expr = minE)
    expressed_genes <- row.names(subset(fData(mono.set), num_cells_expressed >= minN))
    metrics.tab$gene_universe[i]<-length(expressed_genes)

    save(mono.set,file=paste0("minT",minTi,".mono.set.RData"))

    metrics.tab$medGPC[i]<-median(pData(mono.set)$num_genes_expressed,na.rm=TRUE)
 
    print(paste0(i,"_processed"))

    }

save(metrics.tab,file="metrics.tab.RData")

plotdat<-metrics.tab

ggplot(plotdat)+geom_line(aes(x=minT,y=medGPC),size=1)+theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
ggsave("medGPCvsminT.downscaled.png")

ggplot(plotdat)+geom_line(aes(x=minT,y=gene_universe),size=1)+theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
ggsave("gene_universevsminT.downscaled.png")

write.table(metrics.tab,file="metrics.tab.txt",row.names=FALSE,sep="\t",quote=FALSE)



