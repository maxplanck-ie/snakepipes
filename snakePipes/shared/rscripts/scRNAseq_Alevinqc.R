.libPaths(R.home("library"))
library(alevinQC)

indir <- snakemake@params[["indir"]]
outdir <- snakemake@params[["outdir"]]
samID <- snakemake@params[["samid"]]
outfile <- snakemake@params[["outfile"]]

alevinQCReport(baseDir = indir,
               sampleId = samID, 
               outputFile = outfile, 
               outputFormat = "html_document",
               outputDir = outdir, forceOverwrite = TRUE)
