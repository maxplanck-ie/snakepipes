## Usage: cat DESeq2.R | /package/R-3.2.0/bin/R --vanilla --quiet --args setup.tsv counts.txt 0.05 species.gene_names

# args 1 : design matrix (tsv)
# args 2 : counts.txt
# args 3  : FDR
# args 4 : ENSID --> gene symbol mapping file
# args 5 : path to DE_functions
# args 6 : T/F whether or not the workflow is allele-sepecific
# args 7 : tx2gene file for salmon --> DESeq mode

args = commandArgs(TRUE)


## re invented RNaseq workflow
sampleInfoFilePath <- args[1]
countFilePath <- args[2]
fdr <- as.numeric(args[3])
geneNamesFilePath <- args[4]
importfunc <- args[5]
allelic_info <- args[6]
## if output is from salmon then tx2gene file should be present
tx2gene_file <- args[7]
if(file.exists(tx2gene_file)) {
  tximport <- TRUE
} else {
  stop(paste0("file ", tx2gene_file, " doesn't exist!") )
  tximport <- FALSE
}

topN <- 50
## include functions
library(ggplot2)
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
cat(paste("Salmon --> DESeq2 : ", tximport, "\n"))
cat(paste("Salmon tx2gene file : ", tx2gene_file, "\n"))

## ~~~~~ 1. SETUP ~~~~~
## sampleInfo (setup of the experiment)
sampleInfo <- read.table(sampleInfoFilePath, header = TRUE, stringsAsFactor = F)

## add X at the beginning of rows beginning with a number (makes it consistent to column names of of the count matrix!)
if ( any(grepl("^[0-9]", sampleInfo$name)) ) {
	sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name <- paste("X", sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name, sep="")
}

## ~~~~~~ 2. Check if data is in proper order  ~~~~~
if(isTRUE(tximport)) {
  ## Get data from salmon using TXIMPORT 
  ## get gene annotation
  tx2gene <- read.delim(tx2gene_file, header = FALSE)
  tx2gene <- tx2gene[c(1,2)]
  # check setup table and import
  countdata <- checktable(sample_info = sampleInfo, salmon_dir = dirname(countFilePath), tx2gene_annot = tx2gene)
} else {
  ## Get data from featurecounts
  ## counts
  countdata <- read.table(countFilePath, header=TRUE, check.names = FALSE)
  # check setup table and grep colnames
  countdata <- checktable(countdata, sample_info = sampleInfo, alleleSpecific = allelic_info)
}

## ~~~~~~~ 3. run DESeq wrapper ~~~~~~~~
seqout <- DESeq_basic(countdata, coldata = sampleInfo, fdr = fdr, from_salmon = tximport)

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


## ~~~~~~ 4. report on versions used ~~~~~
sink("DESeq2.session_info.txt")
sessionInfo()
sink()
