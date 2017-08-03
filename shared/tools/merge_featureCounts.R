library(tools)

args <- commandArgs(trailingOnly=T)

outfile = args[[1]]
infiles = args[2:length(args)]

# find out if the output is allelic
isallelic <- function(x) {
  return( (ncol(x) == 9) & (sum(grepl("genome[1|2]", colnames(x)[8:9])) == 2) )
}

# get relevant cols
get_df <- function(infile) {
  cat(infile, "\n")
  bname = unlist(strsplit(basename(infile), ".", fixed=T))[[1]]
  df = read.table(infile, header=T)

  if(isallelic(df) == TRUE) {
  print("Counts are allele-specific")
  df <- df[,c(1,7:9)]
<<<<<<< HEAD
  colnames(df)[2:4] <- c(paste0(bname,"_all"), paste0(bname, "_genome", 1:2) )
=======
  colnames(df)[2:4] <- c(bname, paste(bname, "genome", 1:2, sep = "_") )
>>>>>>> fd100fe10ba6c94e7eca03dcf95223c113008485
  } else {
    df <- df[,c(1,7)]
    colnames(df)[2] = bname
  }

  return(df)
}

# Merge the dfs
doMerge <- function(x, y){
  df <- merge(x, y, by= "Geneid", all.x= TRUE, all.y= TRUE)
  return(df)
}

all_df = lapply( infiles, function(x) get_df(x) )
merged_df <- Reduce(doMerge, all_df)
rownames(merged_df) = merged_df$Geneid
merged_df$Geneid = NULL

# output
print("Merged count file")
dim(merged_df)
head(merged_df)

write.table(merged_df, file=outfile, quote=F, sep="\t", col.names=NA)
