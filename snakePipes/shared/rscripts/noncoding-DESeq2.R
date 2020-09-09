#!/usr/bin/env Rscript

.libPaths(R.home("library"))

## re invented RNaseq workflow
sampleInfoFilePath = snakemake@input[['sampleSheet']]
countFiles = snakemake@input[['cnts']]
fdr = as.numeric(snakemake@params[['fdr']])
geneNamesFilePath = snakemake@input[['symbol_file']]
topN <- 50

sink(paste0(snakemake@params["outdir"], "/DESeq2.session_info.txt"))

## include functions
library(ggplot2)
library(rmarkdown)
library(knitcitations)
library(DESeq2)
source(paste0(snakemake@scriptdir, "/DE_functions.R"))

## ~~~~~ 1. SETUP ~~~~~
## sampleInfo (setup of the experiment)
sampleInfo = read.table(sampleInfoFilePath, header = TRUE, stringsAsFactors = F)
sampleInfo$condition = as.factor(sampleInfo$condition)
sampleInfo$condition = relevel(sampleInfo$condition, ref = as.character(sampleInfo$condition[1])) # first sample defines base
rownames(sampleInfo) = sampleInfo$name

## ~~~~~~ 2. Read in data ~~~~~
cnts = lapply(sampleInfo$name, function(x) read.delim(paste0("TEcount/", x, ".cntTable")))
cnts2 = as.matrix(sapply(cnts,function(x) x[,2]))
row.names(cnts2) = cnts[[1]][,1]
colnames(cnts2) = sampleInfo$name
gene_cnts = cnts2[grep(":", rownames(cnts2), invert=TRUE),]  # we assume gene IDs lack :
repeat_cnts = cnts2[grep(":", rownames(cnts2)),]

## ~~~~~~~ 3. Gene-level differences ~~~~~~~~
seqout <- DESeq_basic(gene_cnts, coldata = sampleInfo, fdr = fdr, alleleSpecific = FALSE, from_salmon = FALSE)
DESeq_writeOutput(DEseqout = seqout,
                  fdr = fdr, outprefix = paste0(snakemake@params["outdir"], "/genes"),
                  geneNamesFile = geneNamesFilePath)

## ~~~~~~~ 4. Name-level differences ~~~~~~~~
seqout2 <- DESeq_basic(repeat_cnts, coldata = sampleInfo, fdr = fdr, size_factors = sizeFactors(seqout[['dds']]))
DESeq_writeOutput(DEseqout = seqout2,
                  fdr = fdr, outprefix = paste0(snakemake@params["outdir"], "/repeat_name"),
                  geneNamesFile = "")

## ~~~~~~~ 5. Class-level differences ~~~~~~~~
repgrps = sapply(strsplit(row.names(repeat_cnts), ":"), function(x) x[2])
foo = lapply(split(seq_along(repgrps), repgrps),
       function(ind) colSums(repeat_cnts[ind,,drop=FALSE]))[order(unique(repgrps))]
class_counts = do.call(rbind, foo)
seqout3 <- DESeq_basic(class_counts, coldata = sampleInfo, fdr = fdr, size_factors = sizeFactors(seqout[['dds']]))
DESeq_writeOutput(DEseqout = seqout3,
                  fdr = fdr, outprefix = paste0(snakemake@params["outdir"], "/repeat_class"),
                  geneNamesFile = "")

## ~~~~~~~ 6. Family-level differences ~~~~~~~~
repgrps = sapply(strsplit(row.names(repeat_cnts), ":"), function(x) x[3])
foo = lapply(split(seq_along(repgrps), repgrps),
       function(ind) colSums(repeat_cnts[ind,,drop=FALSE]))[order(unique(repgrps))]
family_counts = do.call(rbind, foo)
seqout4 <- DESeq_basic(family_counts, coldata = sampleInfo, fdr = fdr, size_factors = sizeFactors(seqout[['dds']]))
DESeq_writeOutput(DEseqout = seqout4,
                  fdr = fdr, outprefix = paste0(snakemake@params["outdir"], "/repeat_family"),
                  geneNamesFile = "")

## ~~~~~~~ 7.  Create report ~~~~~~~~~~~~~~
bib <- c(
    knitcitations = citation('knitcitations'),
    DT = citation('DT'),
    ggplot2 = citation('ggplot2'),
    knitr = citation('knitr')[3],
    rmarkdown = citation('rmarkdown')[1],
    pheatmap = citation('pheatmap'),
    RColorBrewer = citation('RColorBrewer'),
    DESeq2 = citation('DESeq2'))

write.bibtex(bib, file = 'citations.bib')
file.copy(paste0(snakemake@scriptdir, "/DESeq2Report.Rmd"), to = 'DESeq2_report_basic.Rmd')

## TODO we need 4 of these...
outprefix = "DEseq_basic"
cite_options(citation_format="text", style="html", cite.style="numeric", hyperlink=TRUE)
render('DESeq2_report_basic.Rmd',
       output_file = paste0(snakemake@params["outdir"], "/DESeq2_report_genes.html"),
       output_format = "html_document",
       clean = TRUE,
       params = list(
           DEseqoutRdata = paste0(snakemake@params["outdir"], "/genes_DESeq.Rdata"),
           ddr.df = paste0(snakemake@params["outdir"], "/genes_DEresults.tsv"),
           countdata = paste0(snakemake@params["outdir"], "/genes_counts_DESeq2.normalized.tsv"),
           coldata = sampleInfo,
           fdr = 0.05,
           heatmap_topN = 20,
           geneNamesFile = geneNamesFilePath))

render('DESeq2_report_basic.Rmd',
       output_file = paste0(snakemake@params["outdir"], "/DESeq2_report_repeat_name.html"),
       output_format = "html_document",
       clean = TRUE,
       params = list(
           DEseqoutRdata = paste0(snakemake@params["outdir"], "/repeat_name_DESeq.Rdata"),
           ddr.df = paste0(snakemake@params["outdir"], "/repeat_name_DEresults.tsv"),
           countdata = paste0(snakemake@params["outdir"], "/repeat_name_counts_DESeq2.normalized.tsv"),
           coldata = sampleInfo,
           fdr = 0.05,
           heatmap_topN = 20,
           geneNamesFile = geneNamesFilePath))

render('DESeq2_report_basic.Rmd',
       output_file = paste0(snakemake@params["outdir"], "/DESeq2_report_repeat_class.html"),
       output_format = "html_document",
       clean = TRUE,
       params = list(
           DEseqoutRdata = paste0(snakemake@params["outdir"], "/repeat_class_DESeq.Rdata"),
           ddr.df = paste0(snakemake@params["outdir"], "/repeat_class_DEresults.tsv"),
           countdata = paste0(snakemake@params["outdir"], "/repeat_class_counts_DESeq2.normalized.tsv"),
           coldata = sampleInfo,
           fdr = 0.05,
           heatmap_topN = 20,
           geneNamesFile = geneNamesFilePath))

render('DESeq2_report_basic.Rmd',
       output_file = paste0(snakemake@params["outdir"], "/DESeq2_report_repeat_family.html"),
       output_format = "html_document",
       clean = TRUE,
       params = list(
           DEseqoutRdata = paste0(snakemake@params["outdir"], "/repeat_family_DESeq.Rdata"),
           ddr.df = paste0(snakemake@params["outdir"], "/repeat_family_DEresults.tsv"),
           countdata = paste0(snakemake@params["outdir"], "/repeat_family_counts_DESeq2.normalized.tsv"),
           coldata = sampleInfo,
           fdr = 0.05,
           heatmap_topN = 20,
           geneNamesFile = geneNamesFilePath))

## ~~~~~~ 8. report on versions used ~~~~~
sessionInfo()
sink()

## ~~~~~~ 9. clean up ~~~~~
file.remove('citations.bib')
