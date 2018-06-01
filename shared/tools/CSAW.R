#!/usr/bin/env Rscript
## ChIPseq differential binding workflow

sampleInfoFilePath <- snakemake@input[["sample_info"]]  #"samplesheet.tab"
fdr <- as.numeric(snakemake@params[["fdr"]])
paired <- as.logical(snakemake@params[["paired"]])
fraglength <- as.numeric(snakemake@params[["fragment_length"]])  # This needs to be figured out somehow
windowSize <- as.numeric(snakemake@params[["window_size"]])
importfunc <- snakemake@params[["importfunc"]]  #"DB_functions.R"
allelic_info <- as.logical(snakemake@params[["allele_info"]])

## create output directory

# include functions
sink("CSAW/CSAW.log", append=TRUE)
source(paste(snakemake@config[["main_dir_path"]], snakemake@params[["importfunc"]], sep="/"))
suppressPackageStartupMessages(library(GenomicRanges))
## fix default FDR significance threshold
if ( is.na(fdr) ) fdr <- 0.05
if (!dir.exists("CSAW")) dir.create("CSAW")

## print the info
cat(paste("Working dir:", getwd(), "\n"))
cat(paste("Sample info CSV:", sampleInfoFilePath, "\n"))
cat(paste("FDR:", fdr, "\n"))
cat(paste("paired-end? :", paired, "\n"))
cat(paste("allele-specific? :", allelic_info, "\n"))

## sampleInfo (setup of the experiment)
sampleInfo <- read.table(sampleInfoFilePath, header = TRUE, colClasses = c("character", "factor"))

## is paired end? : define read params
pe = "none"
if(isTRUE(paired)) {
    pe = "both"
    d = read.delim("deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv")
    fraglength = median(d[,6])
}
pe_param <- csaw::readParam(max.frag = 500, pe = pe)  # Some CSAW functions explode the processor count with >1 core

## Read data
chip_object <- readfiles_chip(sampleInfo = sampleInfo,
                              fragment_length = fraglength,
                              window_size = windowSize,
                              alleleSpecific = allelic_info,
                              pe.param = pe_param)

## merge all peaks from the samples mentioned in sampleinfo to test (exclude "control")
# get files to read from MACS
fnames <- sampleInfo[sampleInfo$condition != "control",]$name

allpeaks <- lapply(fnames, function(x) {
    narrow <- paste0("MACS2/",x,".filtered.BAM_peaks.narrowPeak")
    broad <- paste0("MACS2/",x,".filtered.BAM_peaks.broadPeak")
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
chip_object <- tmmNormalize_chip(chip_object, binsize = 10000, plotfile = "CSAW/TMM_normalizedCounts.pdf")

## get DB regions
print("Performing differential binding")
chip_results <- getDBregions_chip(chip_object, plotfile = "CSAW/DiffBinding_modelfit.pdf")

## write output
print("Writing output")
writeOutput_chip(chip_results, outfile_prefix = "CSAW/DiffBinding", fdrcutoff = fdr)


## save data
print("Saving data")
sink()
save(chip_object, chip_results, file = "CSAW/DiffBinding_analysis.Rdata")

#### SESSION INFO

sink("CSAW/CSAW.session_info.txt")
sessionInfo()
sink()

print("DONE..!")
