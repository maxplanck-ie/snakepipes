#run in R3.4.0

.libPaths(R.home("library"))

#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))
#readRenviron("/home/sikora/.Renviron")
options(stringsAsFactors=FALSE,na.rm=TRUE)

require(monocle)
require(Seurat)

set.seed(314)

mtab<-commandArgs(trailingOnly=TRUE)[2]

load(mtab)

metric<-commandArgs(trailingOnly=TRUE)[3]

##select 'best' threshold
metrics.tab<-metrics.tab[order(metrics.tab$minT,decreasing=TRUE),]
minTi<-metrics.tab$minT[which.max(metrics.tab[,metric])]

load(paste0("minT",minTi,".mono.set.RData"))

disp_table <- dispersionTable(mono.set)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
mono.set <- setOrderingFilter(mono.set, unsup_clustering_genes$gene_id)
plot_ordering_genes(mono.set)
ggsave(paste0("mono.set.",minTi,".disp.estim.png"))

mono.set <- reduceDimension(mono.set, max_components = 2,reduction_method = 'tSNE', verbose = T)
mono.set <- clusterCells(mono.set, num_clusters = NULL)
plot_cell_clusters(mono.set,  1, 2, color="Cluster")
ggsave(paste0("mono.set.",minTi,".tsne.auto.Cluster.png"))
rho.sel<-quantile(pData(mono.set)$rho,0.90)
delta.sel<-quantile(pData(mono.set)$delta,0.99)
plot_rho_delta(mono.set,rho_threshold = rho.sel, delta_threshold = delta.sel)
ggsave(paste0("mono.set.",minTi,".rho_delta.png"))
mono.set <- clusterCells(mono.set, num_clusters = NULL,rho_threshold = rho.sel, delta_threshold = delta.sel,skip_rho_sigma=TRUE)
plot_cell_clusters(mono.set,  1, 2, color="Cluster")
ggsave(paste0("mono.set.",minTi,".tsne.thd.Cluster.png"))

save(mono.set,file=paste0("minT",minTi,".mono.set.RData"))

seuset <- exportCDS(mono.set, 'Seurat')

seuset<-SetIdent(seuset, ident.use = pData(mono.set)$Cluster)
save(seuset,file=paste0("minT",minTi,".seuset.RData"))

table(seuset@ident)

markers <- FindAllMarkers(object = seuset,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)

require(dplyr)
top2 <- markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(top2,paste0("minT",minTi,".Top2markers.txt"),sep="\t",row.names=FALSE,quote=FALSE)
DoHeatmap(object = seuset,genes.use = top2$gene,slim.col.label = TRUE,remove.key = TRUE)
ggsave(paste0("minT",minTi,".Top2markers.heatmap.png"))
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
write.table(top10,paste0("minT",minTi,".Top10markers.txt"),sep="\t",row.names=FALSE,quote=FALSE)
DoHeatmap(object = seuset,genes.use = top10$gene,slim.col.label = TRUE,remove.key = TRUE)
ggsave(paste0("minT",minTi,".Top10markers.heatmap.png"))

seuset <- FindVariableGenes(object = seuset,mean.function = ExpMean,dispersion.function = LogVMR,x.low.cutoff = 0.0125,x.high.cutoff = 3,y.cutoff = 0.5,do.plot=FALSE) 
seuset <- RunPCA(object = seuset,pc.genes = seuset@var.genes,do.print = TRUE,pcs.print = 1:5,genes.print = 5)
seuset <- RunTSNE(object = seuset,dims.use = 1:5,do.fast = TRUE) ##just to initiate the slot correctly
seuset@dr$tsne@cell.embeddings<-t(mono.set@reducedDimA)
colnames(seuset@dr$tsne@cell.embeddings)<-c("tSNE_1","tSNE_2")
rownames(seuset@dr$tsne@cell.embeddings)<-rownames(pData(mono.set))
save(seuset,file=paste0("minT",minTi,".seuset.RData"))
TSNEPlot(object = seuset) #group.by="ident",plot.order
ggsave(paste0("minT",minTi,".seuset.tSNE.png"))

for(i in seq_along(unique(top2$cluster))){
    clu<-unique(top2$cluster)[i]
    VlnPlot(object = seuset, features.plot = top2$gene[top2$cluster %in% clu]) ##
    ggsave(paste0("Top2.clu",clu,".violin.png"),width=12,height=6)
    FeaturePlot(seuset,top2$gene[top2$cluster %in% clu],cols.use = c("lightgrey", "blue"),nCol = 2)
    ggsave(paste0("Top2.clu",clu,".featurePlot.png"),width=12,height=6)

}

sink("sessionInfo.txt")
print(sessionInfo())
sink()
