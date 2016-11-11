## Usage: Rscript library_type.R library_type.tsv PE TopHat2 fr-firststrand htseq-count

args = commandArgs(TRUE)

library_type_tsv = args[[1]] # library_type.tsv
library_type = args[[2]] # PE or SE
source = args[[3]] # e.g. TopHat2
source_value = args[[4]] # e.g. fr-firststrand
target = args[[5]] # HISAT2

df = read.table(library_type_tsv, header=T, stringsAsFactors=F, check.names=F)

cat( df[ df["type"]==library_type & df[source]==source_value, ][[target]], "\n" )
