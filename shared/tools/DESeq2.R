## Usage: cat DESeq2.R | /package/R-3.2.0/bin/R --vanilla --quiet --args setup.tsv counts.txt 0.05 species.gene_names

# args 1 : design matrix (tsv)
# args 2 : counts.txt
# args 3  : FDR
# args 4 : ENSID --> gene symbol mapping file
# args 5 : top N genes to plot
# args 6 : T/F whether or not the workflow is allele-sepecific

args = commandArgs(TRUE)

source("DE_functions.R")

## re invented RNaseq workflow
sampleInfoFilePath <- args[1]
countFilePath <- args[2]
geneNamesFilePath <- args[4]
fdr <- as.numeric(args[3])
topN <- as.numeric(args[5])
allelic_info <- args[6]

## fix default FDR significance threshold
if ( is.na(fdr) ) fdr <- 0.05

## fix default topN genes to plot
if ( is.na(topN) ) topN = 50

## print the info
cat(paste("Working dir:", getwd(), "\n"))
cat(paste("Sample info CSV:", sampleInfoFilePath, "\n"))
cat(paste("Count file:", countFilePath, "\n"))
cat(paste("FDR:", fdr, "\n"))
cat(paste("Gene names:", geneNamesFilePath, "\n"))
cat(paste("Number of top N genes:", topN, "\n"))


## sampleInfo (setup of the experiment)
sampleInfo <- read.table(sampleInfoFilePath, header = TRUE, stringsAsFactor = F)

## add X at the beginning of rows beginning with a number (makes it consistent to column names of of the count matrix!)
if ( any(grepl("^[0-9]", sampleInfo$name)) ) {
	sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name <- paste("X", sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name, sep="")
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


## report on versions used
sink("DESeq2.session_info.txt")
sessionInfo()
sink()
