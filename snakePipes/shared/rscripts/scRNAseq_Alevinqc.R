.libPaths(R.home("library"))
library(alevinQC)
args <- commandArgs(trailingOnly=TRUE) 
indir <- args[1]
outdir <- args[2]
samID <- args[3]
outfile <- args[4]

alevinQCReport(baseDir = indir,
               sampleId = samID, 
               outputFile = outfile, 
               outputFormat = "html_document",
               outputDir = outdir, forceOverwrite = TRUE)
