args<-commandArgs(TRUE)

infile = args[[1]]

bname = unlist(strsplit(basename(infile),"[.]"))[[1]]
dname = dirname(infile)

lines = readLines(infile)
s = lines[grep("## HISTOGRAM", lines)+1:length(lines)]
d<-read.delim(textConnection(s), header=T, sep="\t", strip.white=TRUE)
d

read_orientations = colnames(d)[2:length(colnames(d))]
read_orientations

generate_histogram <- function(read_orientation) {
  x = d$insert_size
  y = as.numeric(unlist(d[read_orientation]))
  png(file.path(dname, paste(bname, ".", read_orientation, ".png", sep="")), width=800, height=600, pointsize = 14)
  barplot(y, names.arg=x, main=read_orientation, col="black", border="black", space=0)
  dev.off()
}

sapply(read_orientations, generate_histogram)
