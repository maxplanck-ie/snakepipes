## Usage: cat DESeq2.R | /package/R-3.2.0/bin/R --vanilla --quiet --args setup.tsv counts.txt 0.05 species.gene_names

# args 1 : design matrix (tsv)
# args 2 : counts.txt
# args 3  : FDR
# args 4 : ENSID --> gene symbol mapping file
# args 5 : top N genes to plot
# args 6 : T/F whether or not the workflow is allele-sepecific

# args 1 : design matrix (tsv)
# args 2 : counts.txt
# args 3  : FDR
# args 4 : BioMart file with ensembl and symbol names
# args 5 : top N genes to plot

args = commandArgs(TRUE)
<<<<<<< HEAD

=======
print(args)
>>>>>>> fd100fe10ba6c94e7eca03dcf95223c113008485

## re invented RNaseq workflow
sampleInfoFilePath <- args[1]
countFilePath <- args[2]
fdr <- as.numeric(args[3])
geneNamesFilePath <- args[4]
importfunc <- args[5]
allelic_info <- args[6]
topN <- 50
## include functions
source(importfunc)

## fix default FDR significance threshold
if ( is.na(fdr) ) fdr <- 0.05

## fix default topN genes to plot
if ( is.na(topN) ) topN <- 50

## print the info
cat(paste("Working dir:", getwd(), "\n"))
cat(paste("Sample info CSV:", sampleInfoFilePath, "\n"))
cat(paste("Count file:", countFilePath, "\n"))
cat(paste("FDR:", fdr, "\n"))
cat(paste("Gene names:", geneNamesFilePath, "\n"))
cat(paste("Number of top N genes:", topN, "\n"))

<<<<<<< HEAD
=======
## For debugging only!!! #######################################################
## setwd("output_dir")
## args = c('setup_table.tsv',
##          'counts.txt',
##          '0.05',
##          'species.gene_names')
################################################################################

print("Running DESeq2 from rna-seq-qc...")
## sampleInfo (setupt of the experiment)
sampleInfo = read.table(sampleInfoFilePath, header=TRUE, stringsAsFactor=F)
## add X at the beginning of rows beginning with a number (makes it consistent to column names of of the count matrix!)
if ( any(grepl("^[0-9]", sampleInfo$name)) ) {
  sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name = paste("X", sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name, sep="")
}
sampleInfo = DataFrame(as.data.frame(unclass(sampleInfo)))
# print sample names
as.character(sampleInfo$name)

## count matrix (e.g. from DESeq or featureCounts)
countdata = read.table(countFilePath, header=TRUE)
countdata = DataFrame(round(countdata))
colnames(countdata)

## 1st check: if names of the setup table are subset of the count matrix column names
if ( ! all( is.element(sort(sampleInfo[,1]), sort(colnames(countdata))) ) ) {
  cat("Error! Count table column names and setup table names do NOT match!\n")
  print(as.character(sampleInfo[,1]))
  print(colnames(countdata))
  quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
}

## extract only the columns specified in the setup table
countdata = countdata[,as.character(sampleInfo[,1])]
head(countdata)
colnames(countdata)

## 2nd: check if the ORDER of sample names matches as well
if ( ! all(as.character(sampleInfo$name) == colnames(countdata)) ) {
  cat("Error! Count table column names and setup table names do NOT match!\n")
  print(as.character(sampleInfo$name))
  print(colnames(countdata))
  quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
}

# create DEseq dds object
dds = DESeqDataSetFromMatrix(
  countData = countdata,
  colData = sampleInfo,
  design = ~ condition)

# print
dds

## reorder conditions by sampleInfo
if ( dds$condition[[1]] != levels(dds$condition)[[1]] ) {
  dds$condition = relevel(dds$condition, as.character(dds$condition[[1]]) )
}

colnames(dds) = sampleInfo$name
head(assay(dds))

################################################################################
## Extra data collecting some measures for every sample (e.g. total
## counts, scaling factors,...)
info <- data.frame(row.names=sampleInfo$name)
################################################################################
## counts per sample (not needed , use library sizes)
apply(assay(dds), 2, sum)
info$total_counts = apply(assay(dds), 2, sum)   # add to info
info

## count table
head(assay(dds))
write.table(assay(dds),"counts.tsv", sep="\t", quote=FALSE, col.names=NA) # save to file

## DE analysis
assign("last.warning", NULL, envir = baseenv())
dds = DESeq(dds)
warnings()
sink("DESeq2.WARNINGS.txt"); warnings(); sink() # save warnings to file

## show size factors used for read count normalisation
sizeFactors(dds)
info$size_factors = sizeFactors(dds)
info

# save normalized counts to file
write.table(counts(dds, normalized=T),"DESeq2.counts_normalized.tsv", sep="\t", quote=FALSE, col.names=NA)


## Expression density plot
toplot = data.frame(counts(dds, normalized=T))
toplot = stack(toplot, select=colnames(toplot))
p = ggplot( toplot, aes(values, colour=ind, alpha=0.5))
p + geom_line(aes(color=ind), stat="density", alpha=0.5) +
  scale_x_log10(name="\nnormalized counts", breaks=c(0.1,1,10,100,1000,10000,100000), limits=c(0.1,100000) ) +
  scale_y_continuous(name="density\n") +
  scale_colour_discrete(name="Samples") +
  geom_vline(xintercept=10, colour="grey", linetype = "dashed") +
  theme_minimal() +
  ggtitle("Density plot\n") +
  theme()
ggsave(file=sprintf("Fig8.Density_plot.sample_read_counts.pdf"), width=7, height=6)


# ## size factors plot
# plotdata = data.frame(name=colnames(dds), size_factor=sizeFactors(dds))
# plotdata
# ggplot(plotdata, aes(x=name)) +
#   geom_bar(aes(weight=size_factor), fill="darkseagreen", width=.7) +
#   theme_bw(base_size=14) +
#   geom_text(aes(y=size_factor, label=sprintf("%.3f",size_factor)), size = 4, hjust = 0.5, vjust = -1) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("") +
#   ylab("Size factor")
# ggsave(file=sprintf("Fig7.Size_factors.pdf"))

## dispersion plot
pdf("Fig1.dispersion_plot.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()
>>>>>>> fd100fe10ba6c94e7eca03dcf95223c113008485

## sampleInfo (setup of the experiment)
sampleInfo <- read.table(sampleInfoFilePath, header = TRUE, stringsAsFactor = F)

<<<<<<< HEAD
## add X at the beginning of rows beginning with a number (makes it consistent to column names of of the count matrix!)
if ( any(grepl("^[0-9]", sampleInfo$name)) ) {
	sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name <- paste("X", sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name, sep="")
=======
plotVolcano <- function(res_obj, data=plot) {
  # Volcano plot
  xlim = c(-4,4)
  ylim = c(0,20)
  cex=c(0.3,0.5)
  plotdata = data.frame(log2FoldChange=res_obj$log2FoldChange, padj=res_obj$padj )
  plotdata = plotdata[!is.na(plotdata),]
  plotdata$cex = cex[[1]]
  plotdata$pch = 16
  plotdata$col = "#525252"
  plotdata$col[plotdata$padj<=fdr] = "#cd0000"

  plotdata$pch[plotdata$log2FoldChange<xlim[[1]]] = 5
  plotdata$cex[plotdata$log2FoldChange<xlim[[1]]] = cex[[2]]
  plotdata$log2FoldChange[plotdata$log2FoldChange<xlim[[1]]] = xlim[[1]]

  plotdata$pch[plotdata$log2FoldChange>xlim[[2]]] = 5
  plotdata$cex[plotdata$log2FoldChange>xlim[[2]]] = cex[[2]]
  plotdata$log2FoldChange[plotdata$log2FoldChange>xlim[[2]]] = xlim[[2]]

  plotdata$pch[-log10(plotdata$padj) > ylim[[2]]] = 2
  plotdata$cex[-log10(plotdata$padj) > ylim[[2]]] = cex[[2]]
  plotdata$padj[-log10(plotdata$padj) > ylim[[2]]] = 10^-ylim[[2]]

  #head(plotdata)
  #dim(plotdata)
  plot(plotdata$log2FoldChange, -log10(plotdata$padj),
       main=sprintf("Volcano plot\n(FDR: %.2f, up: %d, down: %d)",fdr,length(de_up[,1]),length(de_down[,1])),
       xlab="log2-fold change",
       ylab="-log10 q-value",
       xlim=xlim,
       ylim=ylim,
       cex=plotdata$cex, pch=plotdata$pch,
       col=plotdata$col)
  abline(h=-log10(fdr), col=rgb(0,0,1,0.5), lwd=4)
  abline(v=0, col=rgb(1,0,0,0.5), lwd=4)
}

################################################################################

################################################################################
## gene names dict if available
################################################################################

if (file.exists(geneNamesFilePath)) {
  cat(paste("Gene names file found\n"))
  #geneNames = read.csv(geneNamesFilePath, sep="\t", header=F, row.names=1, stringsAsFactors=FALSE)
  geneNames = read.csv(geneNamesFilePath, sep="\t", header=F, stringsAsFactors=FALSE)
  geneNames = geneNames[!duplicated(geneNames[,1]),]
  rownames(geneNames) = geneNames[,1]
  geneNames[,1] = NULL
  head(geneNames)

  if (length( intersect( gsub("\\..*", "", res@rownames), rownames(geneNames) ) ) > 0) {
    cat(paste("Names matching to IDs found\n"))

    ## make a dictionary
    gene_names_dic = geneNames[[1]]
    names(gene_names_dic) = rownames(geneNames)
    ##gene_names_dic["ENSMUSG00000025332"]
  }
}

## generate a dataframe from ids
id_to_gene_name = function(ids) {
  d = data.frame(IDs=gsub("\\..*", "", ids), gene_names=NA)
  d$gene_names = gene_names_dic[ as.character(d$IDs) ]
  head(d)

  # some might be NAs; replace those by original ID
  d[which(is.na(d$gene_names)),]$gene_names = as.character(d[which(is.na(d$gene_names)),]$IDs)
  head(d)
  return(d$gene_name)
}

gene_names_df <- function(obj) {
  df = as.data.frame(obj)
  if (dim(df)[[1]] > 0) {
    if ( exists("gene_names_dic") ) {
      df$gene_names = id_to_gene_name(rownames(df))
    }
  }
  return(df)
>>>>>>> fd100fe10ba6c94e7eca03dcf95223c113008485
}

## counts
countdata <- read.table(countFilePath, header=TRUE)
# check setup table
checktable(alleleSpecific = allelic_info)

## now run DESeq wrapper
seqout <- DESeq_basic(countdata, coldata = sampleInfo, fdr = fdr)

DESeq_downstream(DEseqout = seqout, countdata, sampleInfo,
		     fdr = fdr, outprefix = "DEseq_basic", heatmap_topN = topN,
		     geneNamesFile = geneNamesFilePath)

## Run allele-sepecific DESeq wrapper (if asked for)
if (isTRUE(allelic_info)) {
	seqout_allelic <- DESeq_allelic(countdata, coldata = sampleInfo, fdr = fdr)

	DESeq_downstream(DEseqout = seqout_allelic, countdata, sampleInfo,
			     fdr = fdr, outprefix = "DEseq_allelic", heatmap_topN = topN,
			     geneNamesFile = geneNamesFilePath)
	}


<<<<<<< HEAD
## report on versions used
=======

# attr(res,"filterThreshold")
# plot(attr(res,"filterNumRej"), type="b",
#      ylab="number of rejections",
#      xlab="quantiles of mean of normalized counts")

# plot(res$baseMean+1, -log10(res$padj),
#      log="x",
#      xlab="mean of normalized counts",
#      ylab="-log10 padj",
#      cex=.4, col=rgb(0,0,0,.3))
# abline(h=-log10(fdr), col="red", lwd=1)
#
# plot(metadata(res)$filterNumRej,
#      type="b", ylab="number of rejections",
#      xlab="quantiles of filter")
# lines(metadata(res)$lo.fit, col="red")
# abline(v=metadata(res)$filterTheta)

################################################################################
## rlog transform; for clustering and ordination (e.g PCA)
rld = rlog(dds)
head(assay(rld))

# save rlog tranformed counts to file
write.table(gene_names_df(assay(rld)),"DESeq2.counts_rlog.tsv", sep="\t", quote=FALSE, col.names=NA)

# show that DEseq's rlog works; not really needed for data analysis
##par(mfrow=c(1,2))
##plot(log2(1+counts(dds, normalized=T)[,1:2]), col="#00000020", pch=20, cex=0.3)     #log2
##plot(assay(rld)[,1:2], col="#00000020", pch=20, cex=0.3)                            #rlog (DESeq); is superior

## Sample distances
sampleDists = dist(t(assay(rld)))
sampleDists

## Euclidean sample distance heatmap
sampleDistMatrix = as.matrix(sampleDists)
sampleDistMatrix
rownames(sampleDistMatrix) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
colnames(sampleDistMatrix) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
sampleDistMatrix

colours = colorRampPalette(rev(brewer.pal(9, "GnBu")))(255)
pdf("Fig5.Heatmap.pdf", width=6, height=6)
par(cex.main=1)
heatmap.2(sampleDistMatrix,trace="none",col=colours,
          main="Heatmap\n(Euclidean distances)",
          keysize=1.2,
          notecex=1.1,
          cexRow=1.0, cexCol=1.0, margins=c(10,10),
          cellnote=round(sampleDistMatrix,1),
          notecol="black")
dev.off()

## PCA
data <- plotPCA(rld, intgroup=c("name", "condition"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=name)) +
  geom_hline(aes(yintercept=0), colour="grey") +
  geom_vline(aes(xintercept=0), colour="grey") +
  geom_point(size=5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA\n") +
  scale_shape_manual(values=c(0:18,33:17))
ggsave(file=sprintf("Fig6.PCA.pdf"), width=7, height=6)

##pdf("PCA.pdf")
##plotPCA(rld, intgroup=c("name", "condition"))
##dev.off()

# topN genes by pvalue
if (length(de_total[,1]) > 0) {
  d = data.frame(id=rownames(de_total), padj=de_total$padj)
  if ( length(rownames(d)) < topN ) topN = length(rownames(d))

  d_topx_padj = d[order(d$padj, decreasing=F),][1:topN,]
  d_topx_padj
  plotdata = assay(rld)[as.character(d_topx_padj$id),]  # <- error
  plotdata

  ## test
  setdiff( as.character(d_topx_padj$id), rownames(plotdata))

  # rownames(plotdata) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
  # colnames(plotdata) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")

  if ( exists("gene_names_dic") ) rownames(plotdata) = id_to_gene_name(rownames(plotdata))  # exchange ids by gene names
  plotdata

  pdf(sprintf("Fig7.gene_clustering_top%i_DE_genes.pdf",topN), pointsize = 9)
  heatmap.2(plotdata, scale="row", trace="none", dendrogram="column",
            col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
            main=sprintf("Top %d DE genes\n(by adj. p-value)", topN), keysize=1,
            margins = c(10,10),
            cexRow=0.7, cexCol=0.9)
  dev.off()
}

################################################################################

## cleanup
if ( file.exists("Rplots.pdf") ) { file.remove("Rplots.pdf") }

# report on versions used
>>>>>>> fd100fe10ba6c94e7eca03dcf95223c113008485
sink("DESeq2.session_info.txt")
sessionInfo()
sink()
