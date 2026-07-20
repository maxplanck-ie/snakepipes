#!/usr/bin/env Rscript

.libPaths(R.home("library"))

##instantiate objects from snakemake

merged_peaks<-snakemake@params[["merged_peaks"]]
wdir<-snakemake@params[["outdir"]]

setwd(wdir)

## load libraries
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(GenomicRanges))

## run the code

peaktable<-read.table(merged_peaks)
peaksGR<-GRanges(seqnames=peaktable$V1,ranges=IRanges(start = peaktable$V2,end=peaktable$V3,))
  
peakran<-randomizeFeature(reduce(peaksGR),stranded=TRUE,keep.strand.prop=TRUE,keep.chrom=TRUE,seed=123)
rtracklayer::export(peakran,con=sub("_peak_intersect_f0.5.bed","_randomized_peaks.bed",basename(merged_peaks)))


