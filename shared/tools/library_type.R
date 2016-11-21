## Usage: Rscript library_type.R library_type.tsv PE fr-firststrand TopHat2 htseq-count

args = commandArgs(TRUE)

library_type_tsv = args[[1]] # library_type.tsv
paired = args[[2]] # PE or SE
from_library_type = args[[3]] # e.g. fr-firststrand
from_prg = args[[4]] # e.g. TopHat2
to_prg = args[[5]] # HISAT2

df = read.table(library_type_tsv, header=T, stringsAsFactors=F, check.names=F)

cat( df[ df["type"]==paired & df[from_prg]==from_library_type, ][[to_prg]] )
