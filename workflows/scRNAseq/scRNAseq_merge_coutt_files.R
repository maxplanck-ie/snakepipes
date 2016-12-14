##########################################################################
## script for merging several coutt files with single cell data 
## and output the combined data to file
##
## Author: Steffen Heyne MPI-IE Freiburg

## call example: /package/R-3.3.1/bin/Rscript my_scRNA_path/ beta_cells_28wks.merged.coutt.tsv 1

args <- commandArgs(trailingOnly=TRUE)

## path with all the *.coutt.csv files that should be combined
coutt_path = ""
if (!is.na(args[1])) {
	coutt_path <- args[1] 
} else {
	print("Path param missing! Exit...\n")
	quit(1)
}

## name of merged output file
out_file="merged_coutt.csv"
if (!is.na(args[2])) {
	out_file <- args[2] 
} else {
	print("Output filename param is missing! Exit...\n")
	quit(1)
}

## each coutt uses only 96 from 192 cellbarcodes, so if split is TRUE then identify which set is used and only merge these ones 
split=FALSE
if (!is.na(args[3])) {
	split <- as.logical(as.numeric(args[3])) 
} else {
	print("Split (0/1) is missing! Exit...\n")
	quit(1)
}

######################################################

print(paste(coutt_path,"/*.coutt.csv",sep=""))

dat=data.frame(file=Sys.glob(paste(coutt_path,"/*.coutt.csv",sep="")),stringsAsFactors=F)

## create sample names from file names
dat$name=gsub(".*/(.*).coutt.csv","\\1",dat$file)

print(dat$file)
print(dat$name)

alldata=data.frame(GENEID=NA)

##iterate over all files and combine them into alldata
for (i in 1:length(dat$file)) {
  print(dat$file[i]);
  tmp=read.csv(dat$file[i],sep="\t",header=T)
  rownames(tmp)=tmp$GENEID
  tmp<-tmp[,-1]
  
  ## check if first 96 cells-barcodes or second 96 cells-barcodes were in current sample
  s1=sum(colSums(tmp[,1:96]))
  s2=sum(colSums(tmp[,97:192]))
  if (split & (s1 > s2)) { print("1"); print(s1); tmp<-tmp[,c(seq(1,96))]} else {print("2"); print(s2);tmp<-tmp[,-c(seq(1,96))]}
  
  ## add sample name as prefix to each cell
  colnames(tmp)<-sapply(colnames(tmp), function(x,t) paste0(t,sub("X","_",x)),t=dat$name[i],simplify="array")
  tmp$GENEID<-rownames(tmp)
  alldata<-merge(alldata,tmp,by="GENEID",all=T)
}

## replace NAs with 0
alldata[is.na(alldata)] <- 0
dim(alldata)

## write out table 
write.table(alldata,out_file,sep="\t",col.names=T,quote=F,row.names=F)
