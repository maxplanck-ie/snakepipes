## ChIPseq differential binding workflow

## Usage: cat CSAW.R | /package/R-3.2.0/bin/R --vanilla --quiet --args samplesheet.tsv 0.05 TRUE FALSE

# args 1 : sample information matrix (tsv)
# args 2 : FDR cutoff
# args 3 : paired-end info
# args 4 : seq fragment length
# args 5 : window size to count reads
# args 6 : allele-specific info
# args 7 : path to DB_functions

## get Args
args = commandArgs(TRUE)

sampleInfoFilePath <- args[1] #"samplesheet.tab"
fdr <- as.numeric(args[2])
paired <- as.logical(args[3])
fraglength <- as.numeric(args[4])
windowSize <- as.numeric(args[5])
importfunc <- args[6] #"~/programs/snakemake_workflows/shared/tools/DB_functions.R"
allelic_info <- as.logical(args[7])

# include functions
source(importfunc)
library(GenomicRanges)
## fix default FDR significance threshold
if ( is.na(fdr) ) fdr <- 0.05

## print the info
cat(paste("Working dir:", getwd(), "\n"))
cat(paste("Sample info CSV:", sampleInfoFilePath, "\n"))
cat(paste("FDR:", fdr, "\n"))
cat(paste("paired-end? :", paired, "\n"))
cat(paste("allele-specific? :", allelic_info, "\n"))

## create output directory
dir.create("CSAW")

## sampleInfo (setup of the experiment)
sampleInfo <- read.table(sampleInfoFilePath, header = TRUE, colClasses = c("character", "factor"))

## is paired end? : define read params
if(isTRUE(paired)) {
	pe_param <- csaw::readParam(max.frag = 500, pe = "both")#, restrict = "X")
} else {
	pe_param <- csaw::readParam(max.frag = 500, pe = "none")
}

## Read data
chip_object <- readfiles_chip(sampleInfo = sampleInfo,
					fragment_length = fraglength,
					window_size = windowSize,
					alleleSpecific = allelic_info,
					pe.param = pe_param)

## make QC plot for one sample
first_bam <- head(SummarizedExperiment::colData(chip_object$windowCounts)$bam.files, n = 1)
last_bam <- tail(SummarizedExperiment::colData(chip_object$windowCounts)$bam.files, n = 1)

print(paste0("Making QC plots for first sample : ", first_bam))
makeQCplots_chip(bam.file = first_bam, outplot = "CSAW/QCplots_first_sample.pdf", pe.param = pe_param)

print(paste0("Making QC plots for last sample : ", last_bam))
makeQCplots_chip(bam.file = last_bam, outplot = "CSAW/QCplots_last_sample.pdf", pe.param = pe_param)

## merge all peaks from the samples mentioned in sampleinfo to test
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
		stop("MACS2 output doesn't exist. Nerither ", narrow, " , nor ", broad)
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
save(chip_object, chip_results, file = "CSAW/DiffBinding_analysis.Rdata")

#### SESSION INFO

sink("CSAW/CSAW.session_info.txt")
sessionInfo()
sink()

print("DONE..!")
