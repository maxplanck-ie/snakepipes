.libPaths(R.home("library"))

library(tools)

args <- commandArgs(trailingOnly=T)
# args = c("Name", "TPM",
#          "/data/processing/kilpert/test/snakemake-workflows/RNA-seq/output_PE_100k/Salmon/counts.tsv",
#          "/data/processing/kilpert/test/snakemake-workflows/RNA-seq/output_PE_100k/Salmon/A277_21_WT.quant.sf",
#          "/data/processing/kilpert/test/snakemake-workflows/RNA-seq/output_PE_100k/Salmon/A277_22_WT.quant.sf")

selstring1 = args[[1]]
selstring2 = args[[2]]
outfile = args[[3]]
infiles = args[4:length(args)]

## infile = infiles[[1]]
get_df <- function(infile) {
  cat(infile, "\n")
  bname = unlist(strsplit(basename(infile), ".", fixed=T))[[1]]
  df = read.table(infile, header=T)[c(selstring1, selstring2)] #[,c(1,5)]
  colnames(df)[2] = bname
  return(df)
}

doMerge <- function(x, y){
  df <- merge(x, y, by="Name", all.x= TRUE, all.y= TRUE)
  return(df)
}

all_df = lapply( infiles, function(x) get_df(x) )
merged_df <- Reduce(doMerge, all_df)
rownames(merged_df) = merged_df[,1]
merged_df[,1] = NULL
head(merged_df)
dim(merged_df)

write.table(merged_df, file=outfile, quote=F, sep="\t", col.names=NA)
