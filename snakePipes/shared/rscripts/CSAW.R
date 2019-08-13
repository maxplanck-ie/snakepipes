#!/usr/bin/env Rscript
## ChIPseq differential binding workflow
.libPaths(R.home("library"))

sampleInfoFilePath <- snakemake@input[["sampleSheet"]]  #"samplesheet.tab"
insert_size_metrics <- snakemake@params[["insert_size_metrics"]] # bamPEFragmentSize output
fdr <- as.numeric(snakemake@params[["fdr"]])
lfc <- as.numeric(snakemake@params[["absBestLFC"]])
pairedEnd <- as.logical(snakemake@params[["pairedEnd"]])
fraglength <- as.numeric(snakemake@params[["fragmentLength"]])  # used when the data is not paired end
windowSize <- as.numeric(snakemake@params[["windowSize"]])
importfunc <- snakemake@params[["importfunc"]]  #"DB_functions.R"
allelic_info <- as.logical(snakemake@params[["allele_info"]])
outdir<-snakemake@params[["outdir"]]
yaml_path<-snakemake@params[["yaml_path"]]

##set up a primitive log
logfile <- file(snakemake@log[["err"]], open="wt")
sink(logfile, type="message")


## create output directory

# include functions
#sink("CSAW/CSAW.log", append=TRUE)
source(paste(snakemake@config[["baseDir"]], snakemake@params[["importfunc"]], sep="/"))
suppressPackageStartupMessages(library(GenomicRanges))
## fix default FDR significance threshold
if ( is.na(fdr) ) fdr <- 0.05
if (!dir.exists(outdir)) dir.create(outdir)
setwd(outdir)

## print the info
cat(paste("Working dir:", getwd(), "\n"))
cat(paste("Sample info CSV:", sampleInfoFilePath, "\n"))
cat(paste("FDR:", fdr, "\n"))
cat(paste("LFC:", lfc, "\n"))
cat(paste("paired-end? :", pairedEnd, "\n"))
cat(paste("allele-specific? :", allelic_info, "\n"))

## sampleInfo (setup of the experiment)
sampleInfo <- read.table(sampleInfoFilePath, header = TRUE, colClasses = c("character", "factor"))
## is paired end? : define read params
pe = "none"
if(isTRUE(pairedEnd)) {
    pe = "both"
    d = read.delim(insert_size_metrics)
    fraglength = median(d[,6])
}
pe_param <- csaw::readParam(max.frag = 500, pe = pe)  # Some CSAW functions explode the processor count with >1 core

## Read data
##filter out input using yaml
library(yaml)
y<-read_yaml(yaml_path)
input_list<-unique(unlist(lapply(y[[1]],function(X)X[["control"]])))
if(!is.null(input_list)&&!(input_list=="")){
    sampleInfo<-subset(sampleInfo,!(name %in% input_list))
}


chip_object <- readfiles_chip(sampleSheet = sampleInfo,
                              fragmentLength = fraglength,
                              window_size = windowSize,
                              alleleSpecific = allelic_info,
                              pe.param = pe_param)

## make QC plot for first and last sample
if(pairedEnd=="True"){
    first_bam <- head(SummarizedExperiment::colData(chip_object$windowCounts)$bam.files, n = 1)
    last_bam <- tail(SummarizedExperiment::colData(chip_object$windowCounts)$bam.files, n = 1)

    print(paste0("Making QC plots for first sample : ", first_bam))
    makeQCplots_chip(bam.file = first_bam, outplot = "QCplots_first_sample.pdf", pe.param = pe_param)

    print(paste0("Making QC plots for last sample : ", last_bam))
    makeQCplots_chip(bam.file = last_bam, outplot = "QCplots_last_sample.pdf", pe.param = pe_param)}else{message("No QC plots available for single end reads.")}

## merge all peaks from the samples mentioned in sampleinfo to test (exclude those with 'False' in the UseRegions column)
# get files to read from MACS
if (!is.null(sampleInfo$UseRegions)){
    fnames <- sampleInfo$name[as.logical(sampleInfo$UseRegions)]} else {fnames<-sampleInfo$name}

allpeaks <- lapply(fnames, function(x) {
    narrow <- paste0("../MACS2/",x,".filtered.BAM_peaks.narrowPeak")
    broad <- paste0("../MACS2/",x,".filtered.BAM_peaks.broadPeak")
    # first look for narrowpeak then braod peak
    if(file.exists(narrow)) {
        bed <- read.delim(narrow, header = FALSE)
    } else if (file.exists(broad)) {
        bed <- read.delim(broad, header = FALSE)
    } else {
        stop("MACS2 output doesn't exist. Neither ", narrow, " , nor ", broad)
    }

    bed.gr <- GRanges(seqnames = bed$V1, ranges = IRanges(start = bed$V2, end = bed$V3), name = bed$V4)
    return(bed.gr)
    })
# merge
allpeaks <- Reduce(function(x,y) GenomicRanges::union(x,y), allpeaks)

## keep only these peaks for testing DB
print(paste0("Filtering windows using MACS2 output : ", length(allpeaks) , " regions used (Union of peaks)"))

keep <- overlapsAny(SummarizedExperiment::rowRanges(chip_object$windowCounts), allpeaks)
chip_object$windowCounts <- chip_object$windowCounts[keep,]

## TMM normalize
print("Normalizing using TMM (using 10kb background counts)")
chip_object <- tmmNormalize_chip(chip_object, binsize = 10000, plotfile = "TMM_normalizedCounts.pdf")

## get DB regions
print("Performing differential binding")
chip_results <- getDBregions_chip(chip_object, plotfile = "DiffBinding_modelfit.pdf")

## write output
print("Writing output")
writeOutput_chip(chip_results, outfile_prefix = "DiffBinding", fdrcutoff = fdr,lfccutoff=lfc)


## save data
print("Saving data")
#sink()
save(chip_object, chip_results, file = "DiffBinding_analysis.Rdata")

sink(type="message")
close(logfile)
#### SESSION INFO

sink("CSAW.session_info.txt")
sessionInfo()
sink()

print("DONE..!")
