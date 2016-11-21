library(tools)

args <- commandArgs(trailingOnly=T)

outfile = args[[1]]
infiles = args[2:length(args)]


get_df <- function(infile) {
  cat(infile, "\n")
  bname = unlist(strsplit(basename(infile), ".", fixed=T))[[1]]
  df = read.table(infile, header=T)[,c(1,7)]
  colnames(df)[2] = bname
  return(df)
}

doMerge <- function(x, y){
  df <- merge(x, y, by= "Geneid", all.x= TRUE, all.y= TRUE)
  return(df)
}

all_df = lapply( infiles, function(x) get_df(x) )
merged_df <- Reduce(doMerge, all_df)
rownames(merged_df) = merged_df$Geneid
merged_df$Geneid = NULL
head(merged_df)
dim(merged_df)

write.table(merged_df, file=outfile, quote=F, sep="\t", col.names=NA)
