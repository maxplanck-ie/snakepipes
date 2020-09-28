#this is a modified copy of https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/sce_helpers.R and https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/summarize_spliced_unspliced_dentate_gyrus.R

sink(snakemake@log[["out"]])
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
  cdna_introns <- tximeta(coldata = data.frame(
    names = sampleid,
    files = file.path(alevindir,s,"alevin", "quants_mat.gz"),
    stringsAsFactors = FALSE
  ), type = "alevin")
  uidx <- grep("\\.*I\\.*$", cdna_introns[,2])
  sidx <- grep("\\.*I\\.*$", cdna_introns[,2], invert = TRUE)
  ucounts <- assay(cdna_introns, "counts")[uidx, ]
  scounts <- assay(cdna_introns, "counts")[sidx, ]
  rownames(ucounts) <- gsub("\\.*I\\.*$", "", rownames(ucounts))
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
alevindir<-snakemake@params[["alevindir"]]
samplenames<-snakemake@params[["samplenames"]]
outfile<-snakemake@params[["outfile"]]

## cDNA/introns quantified jointly
tx2gene <- read.table(t2g,header=FALSE,sep="\t",quote="",as.is=TRUE)

sce<- do.call(cbind, lapply(samplenames, function(s) {
        tmp <- read_alevin_cdna_introns(
          alevindir = file.path(alevindir),
          sampleid = s, tx2gene = tx2gene
        )
        colnames(tmp) <- paste0(s, "__", colnames(tmp))
        tmp
      }))

saveRDS(sce,outfile)  

message('done all')
sink()

sink("sessionInfo.txt")
sessionInfo()
sink()


