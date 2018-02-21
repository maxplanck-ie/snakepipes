##########################################################################
## script for merging several coutt files with single cell data 
## and output the combined data to file
##
## Author: Steffen Heyne MPI-IE Freiburg

## call example: /package/R-3.3.1/bin/Rscript my_scRNA_path/ beta_cells_28wks.merged.coutt.tsv 1

require(gtools)
require(dplyr)

args <- commandArgs(trailingOnly=TRUE)

## path with all the *.coutt.csv files that should be combined
coutt_path = NULL
if (!is.na(args[1])) {
  coutt_path <- args[1] 
} else {
  print("Path param missing! Exit...\n")
  stop(1)
}

## name of merged output file
#out_file="merged_coutt.csv"
out_file=NULL
if (!is.na(args[2])) {
	out_file <- args[2] 
} else {
	print("Output filename param is missing! Exit...\n")
	stop(1)
}

## name of file with final cell names 
#cell_names_out_file = "merged_coutt.cell_names.csv"
cell_names_out_file = NULL
if (!is.na(args[3])) {
	cell_names_out_file <- args[3] 
} else {
	print("Cell names output filename param is missing! Exit...\n")
	stop(1)
}

## Are only 96 cells instead of 192 are used per samples? 
## This is ignored if cell_names_path is provided!
is_split_library = FALSE

if (!is.na(args[4])) {
	is_split_library <- as.logical(args[4]) 
} else {
	print("Is split library boolean param is missing! Exit...\n")
	stop(1)
}

## path to cell_names file 
cell_names_path = NULL
if (!is.na(args[5])) {
  cell_names_path <- args[5] 
} else {
  print("Cell names file param missing!\n Use all cells from found coutt files!")
#  quit(1)
}

######################################################

#cell_names_path = "A732_cell_names.tsv"
print(paste(coutt_path,"/*.coutt.csv",sep=""))

dat=data.frame(file=mixedsort(Sys.glob(paste(coutt_path,"/*.coutt.csv",sep=""))),stringsAsFactors=F)
#dat=data.frame(file=mixedsort(Sys.glob(paste("A732_mapping.gencode.genome/Counts/","/Sebastian1*.coutt.csv",sep=""))),stringsAsFactors=F)

## create sample names from file names
dat$name=gsub(".*/(.*).coutt.csv","\\1",dat$file)

print(dat$file)
print(dat$name)

cell_names <- NULL
if (is.null(cell_names_path)){
	print("use all cells in found coutt files!")
} else if (file.exists(cell_names_path) & file.info(cell_names_path)$isdir == FALSE) {
	## cell_names_path is file with ranges of cell names (column names: cell_idx_start,cell_idx_end,cell_prefix)
	print(paste("Read in csv file with cell_names as ranges (",cell_names_path,")"))
	cell_names = read.table(cell_names_path, sep = "\t", stringsAsFactors = F, header = T)
} else if ( file.exists(cell_names_path) & file.info(cell_names_path)$isdir == TRUE ) {
	## cell_names_path is given as directory, use all *.csv files for cell names
	print(paste("Read in all csv files with cell_names from path",cell_names_path))
	files <- list.files(cell_names_path, pattern=".csv$", full.names = TRUE)
	
	if (length(files) == 0){
		print(paste("no csv files found under",cell_names_path))
		stop(1)
	}
	
	for (f in files) {
		dat_csv <- read.csv(f, header=T, sep="\t", stringsAsFactors = F)
		
		if (!all(c("sample","cell_name") %in% names(dat_csv))) {
			print("required columns not found!")
			stop(1)
		}
	#	str(dat_csv)
		cell_names <- rbind(cell_names, dat_csv)
		
	}
} else {
	print(paste("Provided path with cell names does not exists! Exit! (",cell_names_path,")"))
	stop(1)
}


if (all(c("cell_prefix","cell_idx_start") %in% names(cell_names))) {
	cell_names <- cell_names %>% 
		group_by(cell_prefix,cell_idx_start,cell_idx_end) %>% 
		slice(rep(1:n(), each=cell_idx_end-cell_idx_start +1)) %>% 
		mutate(cell_idx=cell_idx_start+row_number(cell_prefix)-1, cell_name = paste(cell_prefix,"_",cell_idx,sep="")) %>% 
		ungroup() %>% 
		arrange(plate,library,cell_idx_start) %>% 
		as.data.frame()
}

str(dat)
alldata=data.frame(GENEID=NA)
## defaults in case no cell names are provided, ie. all cell from each sample/library are used
cell_names_tmp = NULL
libs_per_plate = 2
if (is_split_library) libs_per_plate = 4

for (i in 1:length(dat$file)) {
	
	print(dat$file[i]);
	
	## get all cell names for current sample/library/coutt file
	## empty if no cells (sample not present in cell_names_path) or all cells (cell_names_path = NULL) wanted
	subset <- cell_names[cell_names$sample==dat[i,"name"],]
	if 	(is.data.frame(subset) && nrow(subset)==0) { subset <- NULL }
	if ( is.null(subset) && !is.null(cell_names_path) ){
		print(paste(dat[i,"name"]," ignored, nothing found in cell_names file for that sample!"))
		next
	}
	#str(cell_names)

  tmp <- read.csv(dat$file[i],sep="\t",header=T)
  
  ## rename all cells 
  if (is.null(subset)){
  	colnames(tmp)[2:ncol(tmp)]<-sapply(colnames(tmp)[2:ncol(tmp)], function(x,t) paste0(t,sub("X","_",x)),t=dat$name[i],simplify="array")
  	cell_idx_offset = 0;
	if (is_split_library){
  		s1=sum(colSums(tmp[,2:97]))
  		s2=sum(colSums(tmp[,98:193]))
  		if (s1 > s2) { 
			print("use 1st half"); print(s1); print(s2); 
			tmp<-tmp[,c(1,seq(2,97))]; 
			cell_idx_offset = 0 } 
		else {
			print(" use 2nd half"); print(s2); print(s1); 
			tmp<-tmp[,-c(seq(2,97))];
			cell_idx_offset = 96
		}
	}
  	
	cell_names_tmp <- rbind(cell_names_tmp,data.frame(sample=dat$name[i],plate=((i-1) %/% libs_per_plate)+1, library = ((i-1) %% libs_per_plate)+1, cell_idx=seq(1+cell_idx_offset,ncol(tmp)-1+cell_idx_offset),cell_name=colnames(tmp)[2:ncol(tmp)]))
 
  } else {
	print(sum(colSums(tmp[,c((min(subset$cell_idx)+1):(max(subset$cell_idx)+1))])))
	print(sum(colSums(tmp[,-c(1, (min(subset$cell_idx)+1):(max(subset$cell_idx)+1))])))
	rownames(subset) <- subset$cell_idx
  	colnames(tmp)[(min(subset$cell_idx)+1):(max(subset$cell_idx)+1)] <- sapply(colnames(tmp)[(min(subset$cell_idx)+1):(max(subset$cell_idx)+1)], function(x,t) paste0(subset[sub("X","",x),"cell_name"],collapse = "_"),t=subset,simplify="array")
  	tmp <- tmp[,c(1,(min(subset$cell_idx)+1):(max(subset$cell_idx)+1))]
	
	cell_names_tmp <- rbind(cell_names_tmp,subset)
}
  
  alldata<-merge(alldata,tmp,by="GENEID",all=T)
}

## replace NAs with 0
alldata <- alldata[!is.na(alldata$GENEID),]
alldata[is.na(alldata)] <- 0.0

dim(alldata)

## write out table
write.table(alldata,out_file,sep="\t",col.names=T,quote=F,row.names=F)

#if (is.null(cell_names_path)) {
	cell_names <- cell_names_tmp
#}

write.table(cell_names[,c("sample","plate","library","cell_idx","cell_name")],cell_names_out_file,sep="\t",col.names=T,quote=F,row.names=F)

############## old
#write.table(alldata,"../drougard/A558_mapping.GRCm38_gencode.transcriptome.merged_coutt.v2.tsv",sep="\t",col.names=T,quote=F,row.names=F)

#id_cols <- names(sample_names)[grepl(pattern="plate|library|cell_idx|cell_name",x=names(sample_names))]


#colnames(sample_names)<-c("name","plate","start_idx","end_idx")
#colnames(sample_names)[is.na(colnames(sample_names))] <- paste("test_",seq(1,length(colnames(sample_names)[is.na(colnames(sample_names))])),sep="")

#lib_bounds <- sample_names %>% group_by(name) %>% select(name,start_idx,end_idx) %>% summarise(min=min(start_idx),max=max(end_idx))  %>% as.data.frame()

#vars=setdiff(names(sample_names),"cell_idx")

#name_bounds <- sample_names %>% group_by_(.dots=vars) %>% select(sample,plate,library,Cond1,Cond2,Cond3,cell_idx) %>% summarise(min=min(cell_idx),max=max(cell_idx)) %>% as.data.frame()


#sample_names$id <- apply(sample_names[,id_cols],1,paste,collapse="_")
#sample_names %>% group_by(id,rlid=data.table::rleid(id)) %>% select(id,rlid,cell_idx) %>% summarise(min=min(cell_idx),max=max(cell_idx)) %>% arrange(rlid) %>% as.data.frame






#print(paste(coutt_path,"/*.coutt.csv",sep=""))
#
#
#
#
#dat=data.frame(file=Sys.glob(paste(coutt_path,"/*.coutt.csv",sep="")),stringsAsFactors=F)
#
### create sample names from file names
#dat$name=gsub(".*/(.*).coutt.csv","\\1",dat$file)
#
#print(dat$file)
#print(dat$name)
#
#sample_names <- read.table(sample_name_path,sep = "\t")
#
#colnames(sample_names)<-c("name","plate","start_idx","end_idx")
#colnames(sample_names)[is.na(colnames(sample_names))] <- paste("test_",seq(1,length(colnames(sample_names)[is.na(colnames(sample_names))])),sep="")
#
#lib_bounds <- sample_names %>% group_by(name) %>% select(name,start_idx,end_idx) %>% summarise(min=min(start_idx),max=max(end_idx))  %>% as.data.frame()
#
#dat <- merge(dat,lib_bounds,by = "name")
#
#alldata=data.frame(GENEID=NA)
#
#for (i in 1:length(dat$file)) {
#	
#	print(dat$file[i]);
#	
#	subset <- sample_names[sample_names$name==dat[i,"name"],]
#	
#	if (dim(subset)[1]==0) {
#		print("Sample name is missing for file in cell name table!")
#		print(paste("file: ",dat$file[i]))
#		print(paste("missing name: ",dat$name[i]))
#		quit(1)
#	}
#	
#	tmp <- read.csv(dat$file[i],sep="\t",header=T)
#
#	for (j in 1:nrow(subset)) {	
#	
#		print("subset");
#		prefix_subset = dat$name[i]
#		if (ncol(subset)>4){
#			prefix_subset = paste(subset[j,5:ncol(subset)],collapse="_") 
#		}
#	sub_tmp <- tmp[,subset[j,(subset[j,start_idx+1]):(subset[j,(subset[j,end_idx+1])]]
#	colnames(sub_tmp) <- sapply(colnames(sub_tmp), function(x,t) paste0(t,sub("X","_",x)),t=prefix_subset,simplify="array")
#	}
#}
#
###iterate over all files and combine them into alldata
#for (i in 1:length(dat$file)) {
#  print(dat$file[i]);
#  tmp=read.csv(dat$file[i],sep="\t",header=T)
#  rownames(tmp)=tmp$GENEID
#  tmp<-tmp[,-1]
#  
#  ## check if first 96 cells-barcodes or second 96 cells-barcodes were in current sample
#  s1=sum(colSums(tmp[,1:96]))
#  s2=sum(colSums(tmp[,97:192]))
#  if (split & (s1 > s2)) {
#	  print("1"); 
#	  print(s1); 
#	  print(s2); 
#	  tmp<-tmp[,c(seq(1,96))]
#  } else {
#	  print("2");
#	  print(s2);
#	  print(s1);
#	  tmp<-tmp[,-c(seq(1,96))]
#  }
#  
#  ## add sample name as prefix to each cell
#  colnames(tmp)<-sapply(colnames(tmp), function(x,t) paste0(t,sub("X","_",x)),t=dat$name[i],simplify="array")
#  tmp$GENEID<-rownames(tmp)
#  alldata<-merge(alldata,tmp,by="GENEID",all=T)
#}
#
### replace NAs with 0
#alldata[is.na(alldata)] <- 0
#dim(alldata)
#
### write out table 
#write.table(alldata,out_file,sep="\t",col.names=T,quote=F,row.names=F)
