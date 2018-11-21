#run in R3.4.0
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))
#readRenviron("/home/sikora/.Renviron")
options(stringsAsFactors=FALSE,na.rm=TRUE,rgl.useNULL = TRUE)

library(RaceID)
library(ggplot2)
library(reshape2)

set.seed(314)

mtab<-commandArgs(trailingOnly=TRUE)[2]

load(mtab)

metric<-commandArgs(trailingOnly=TRUE)[3]

##select 'best' threshold
metrics.tab<-metrics.tab[order(metrics.tab$minT,decreasing=TRUE),]
minTi<-metrics.tab$minT[which.max(metrics.tab[,metric])]

load(paste0("sc.minT",minTi,".RData"))

sc<-compdist(sc,metric="logpearson",FSelect=TRUE)
sc<-clustexp(sc,rseed=314,FUNcluster="kmedoids")

png(paste0("sc.minT",minTi,".jaccard.png"))
plotjaccard(sc)
dev.off()

sc<-findoutliers(sc)
sc<-comptsne(sc)

save(sc,file=paste0("sc.minT",minTi,".RData"))

png(paste0("sc.minT",minTi,".MeanVar.png"))
plotbackground(sc)
dev.off()

png(paste0("sc.minT",minTi,".tsne.clu.png"))
plotmap(sc,final=FALSE)
dev.off()

res10L<-lapply(unique(sc@cpart),function(X){
        dg<-clustdiffgenes(sc,X,pvalue=.01)
        dgsub<-dg[dg$fc>=2,]
        if(nrow(dgsub)>0){
        dg<-dgsub
        dg<-head(dg,n=10)
        dg$Cluster<-X
        dg$Gene<-rownames(dg)}else{dg<-NULL}
        return(dg)})
top10<-as.data.frame(do.call(rbind,res10L))
top10<-top10[top10$padj<0.05,]
top10<-top10[with(top10, order(Cluster, padj)),]

write.table(top10,paste0("minT",minTi,".Top10markers.txt"),sep="\t",row.names=TRUE,quote=FALSE)

genes <- unique(top10$Gene)
png(paste0("sc.minT",minTi,".Top10markers.heatmap.png"))
plotmarkergenes(sc,genes)
dev.off()

res2L<-lapply(unique(top10$Cluster),function(X){
        head(top10[top10$Cluster %in% X,],n=2)})
top2<-as.data.frame(do.call(rbind,res2L))
top2<-top2[with(top2, order(Cluster, padj)),]
write.table(top2,paste0("minT",minTi,".Top2markers.txt"),sep="\t",row.names=TRUE,quote=FALSE)

genes <- unique(top2$Gene)
png(paste0("sc.minT",minTi,".Top2markers.heatmap.png"))
plotmarkergenes(sc,genes)
dev.off()

for(i in seq_along(unique(top2$Cluster))){
    clu<-unique(top2$Cluster)[i]
    subdat<-as.data.frame(as.matrix(sc@ndata[top2$Gene[top2$Cluster %in% clu],,drop=FALSE])*min(sc@counts)+.1)
    subdat$GeneID<-rownames(subdat)
    plotdat<-melt(subdat,id.vars="GeneID",value.name="NormExpr",variable.name="CellID")
    plotdat$Cluster<-sc@cpart[match(plotdat$CellID,names(sc@cpart))]
    plotdat$Cluster<-factor(plotdat$Cluster,levels=as.character(unique(plotdat$Cluster))[order(as.numeric(unique(plotdat$Cluster)))])
    ggplot(data=plotdat)+geom_violin(aes(x=Cluster,y=NormExpr,fill=Cluster))+geom_boxplot(aes(x=Cluster,y=NormExpr),width=0.1)+ggtitle(paste0("Cluster ",clu))+facet_wrap(~GeneID)
    ggsave(paste0("Top2.clu",clu,".violin.png"),width=12,height=6)
    pdf(paste0("Top2.clu",clu,".featurePlot.pdf"),bg="white",onefile=TRUE)
    plotexpmap(sc,g=rownames(subdat)[1],n=rownames(subdat)[1],logsc=TRUE,fr=FALSE)
    if(nrow(subdat)>1){
    plotexpmap(sc,g=rownames(subdat)[2],n=rownames(subdat)[2],logsc=TRUE,fr=FALSE)}
    dev.off()
    
}

sink("sessionInfo.txt")
print(sessionInfo())
sink()


