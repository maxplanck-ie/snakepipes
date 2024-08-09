#this is a modified copy of https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/sce_helpers.R and https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/summarize_spliced_unspliced_dentate_gyrus.R
.libPaths(R.home("library"))

wdir<-snakemake@params[["wdir"]]
if (!dir.exists(wdir)) dir.create(wdir)
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

suppressPackageStartupMessages({
  library(tximeta)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(BiocSingular)
})



sce_from_scounts_ucounts <- function(scounts, ucounts) {
  ss <- sum(scounts)
  su <- sum(ucounts)
  
  allgenes <- union(rownames(scounts), rownames(ucounts))
  allcells <- union(colnames(scounts), colnames(ucounts))

  scounts <- as.matrix(scounts)
  scounts <- scounts[match(allgenes, rownames(scounts)), 
                     match(allcells, colnames(scounts))]
  scounts[is.na(scounts)] <- 0
  rownames(scounts) <- allgenes
  colnames(scounts) <- allcells
    
  ucounts <- as.matrix(ucounts)
  ucounts <- ucounts[match(allgenes, rownames(ucounts)), 
                     match(allcells, colnames(ucounts))]
  ucounts[is.na(ucounts)] <- 0
  rownames(ucounts) <- allgenes
  colnames(ucounts) <- allcells
  
  stopifnot(all(rownames(ucounts) == rownames(scounts)))
  stopifnot(all(colnames(ucounts) == colnames(scounts)))
  stopifnot(sum(scounts) == ss)
  stopifnot(sum(ucounts) == su)
  
  SingleCellExperiment(
    assays = list(counts = as(scounts, "dgCMatrix"),
                  spliced = as(scounts, "dgCMatrix"),
                  unspliced = as(ucounts, "dgCMatrix"))
  )
}


read_alevin_cdna_introns <- function(alevindir, sampleid, tx2gene) {
  cdna_introns <- tximport::tximport(
    files = file.path(alevindir,sampleid,"alevin", "quants_mat.gz"),
    type = "alevin",dropInfReps=TRUE,tx2gene=tx2gene)[["counts"]]
  uidx <- grep("\\.*I\\.*$", rownames(cdna_introns))
  sidx <- grep("\\.*I\\.*$", rownames(cdna_introns), invert = TRUE)
  ucounts <- cdna_introns[uidx, ]
  scounts <- cdna_introns[sidx, ]
  rownames(ucounts) <- gsub("\\.*-I\\.*$", "", rownames(ucounts))
  rownames(scounts) <- gsub("\\.*$", "", rownames(scounts))
  cdna_introns <- sce_from_scounts_ucounts(scounts, ucounts)
  rownames(cdna_introns) <- scater::uniquifyFeatureNames(
    ID = rownames(cdna_introns),
    names = tx2gene$gene_name[match(rownames(cdna_introns), tx2gene$gene_id)]
  )
  cdna_introns
}


###snakameke params here
t2g<-snakemake@params[["t2g"]]
g2s<-snakemake@params[["g2s"]]
alevindir<-snakemake@params[["alevindir"]]
samplenames<-snakemake@params[["samplenames"]]
outfile<-snakemake@params[["outfile"]]

print(alevindir)
print(samplenames)

## cDNA/introns quantified jointly
tx2gene <- read.table(t2g,header=FALSE,sep="\t",quote="",as.is=TRUE)
colnames(tx2gene)<-c("transcript_id","gene_id")
gene2symbol<- read.table(g2s,header=FALSE,sep="\t",quote="",as.is=TRUE)
tx2gene$gene_name<-gene2symbol$V2[match(tx2gene$gene_id,gene2symbol$V1)]

sce<- do.call(cbind, lapply(samplenames, function(s) { 
        tmp <- read_alevin_cdna_introns(alevindir = alevindir,sampleid = s, tx2gene = tx2gene)
        colnames(tmp) <- paste0(s, "__", colnames(tmp))
        tmp
      }))

saveRDS(sce,outfile)  

sink("sessionInfo.txt")
sessionInfo()
sink()


