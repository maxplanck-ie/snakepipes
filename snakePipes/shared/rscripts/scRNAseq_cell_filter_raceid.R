#run in R3.4.0
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))


options(stringsAsFactors=FALSE,na.rm=TRUE,rgl.useNULL = TRUE)

library(data.table)
library(ggplot2)
library(RaceID)

#####read-in de-multiplexed count table
csvpath<-commandArgs(trailingOnly=TRUE)[2]

csv.dat<-fread(csvpath,header=TRUE,sep="\t")
#impute 0 where NA
csv.dat[is.na(csv.dat)]<-0

x<-csv.dat
if(!"GENEID" %in% colnames(csv.dat)) { colnames(x)[1]<-"GENEID"} ##for debugging
rownames(x) <- x$GENEID
setkey(x,GENEID)
expdat <- as.matrix(x[grep("ERCC",rownames(x),invert=TRUE),-1,with=FALSE])
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
    
    sc<-SCseq(expdat)
    sc<-filterdata(sc,mintotal=minTi,minexpr=minE,minnumber=minN)
    save(sc,file=paste0("sc.minT",minTi,".RData"))
    
    metrics.tab$num_cells[i]<-ncol(sc@ndata)
    metrics.tab$gene_universe[i]<-length(sc@genes)
    metrics.tab$medGPC[i]<-median(apply(sc@ndata[sc@genes,],2,function(X)sum(X>=minE,na.rm=TRUE)))
    


    print(paste0(i,"_processed"))

    }

save(metrics.tab,file="metrics.tab.RData")

plotdat<-metrics.tab

ggplot(plotdat)+geom_line(aes(x=minT,y=medGPC),size=1)+theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
ggsave("medGPCvsminT.downscaled.png")

ggplot(plotdat)+geom_line(aes(x=minT,y=gene_universe),size=1)+theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
ggsave("gene_universevsminT.downscaled.png")

write.table(metrics.tab,file="metrics.tab.txt",row.names=FALSE,sep="\t",quote=FALSE)



