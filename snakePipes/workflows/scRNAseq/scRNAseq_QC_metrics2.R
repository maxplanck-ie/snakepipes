##########################################################################
## script for merging several coutt files with single cell data 
## and output the combined data to file
##
## Author: Steffen Heyne MPI-IE Freiburg

## call example: /package/R-3.3.1/bin/Rscript my_scRNA_path/ beta_cells_28wks.merged.coutt.tsv 1

require(ggplot2)
require(dplyr)
require(scales)
require(gtools)
require(reshape2)

args <- commandArgs(trailingOnly=TRUE)

## path with all the *.coutt.csv files that should be combined
cellsum_path <- NULL
if (!is.na(args[1])) {
	cellsum_path <- args[1] 
} else {
	print("Cellsum path param missing! Exit...\n")
	quit(1)
}

if (!is.na(args[2])) {
	out_prefix <- args[2] 
} else {
	print("output prefix param missing! Exit...\n")
	quit(1)
}

## Are only 96 cells instead of 192 are used per samples? 
## This is ignored if cell_names_path is provided!
is_split_library = FALSE

if (!is.na(args[3])) {
	is_split_library <- as.logical(args[3]) 
} else {
	print("Is split library boolean param is missing! Exit...\n")
	stop(1)
}

## path to cell_names file 
cell_names_path = NULL
if (!is.na(args[4])) {
	cell_names_path <- args[4] 
} else {
	print("Cell names file param missing!\n Use all cells from found coutt files!")
#  quit(1)
}

## plot format 
plot_format = NULL
if (!is.na(args[5])) {
	plot_format <- args[5] 
} else {
	print("Please provide plot format (pdf,png,None)!")
#  quit(1)
}


###############################################################################

files <- list.files(cellsum_path, pattern=".cellsum",full.names=T)

if (length(files) == 0){
	print(paste("no csv files found under",cellsum_path))
	stop(1)
}

print(files)
sc_dat <- NULL

for (f in files) {
	dat <- read.csv(f, header=T, sep="\t", stringsAsFactors = F)
	
	if (!all(c("sample","cell_idx","READS_UNIQFEAT","UMI") %in% names(dat))) {
		print("required columns not found!")
		stop(1)
	}
	
	sc_dat <- rbind(sc_dat, dat)
}
print(str(sc_dat))

files <- list.files(cellsum_path, pattern=".libsum",full.names=T)

lib_stats <- NULL
if (length(files) != 0){
	
	print(files)
	
	for (f in files) {
		dat <- read.csv(f, header=F, sep="\t", stringsAsFactors = F)
		
		lib_stats <- rbind(lib_stats, dat)
	}

	libstats_reads <- dcast(lib_stats,V1~V2,value.var = "V3")
	libstats_pct <- dcast(lib_stats,V1~V2,value.var = "V4")
	
	#print(str(libstats_reads))
	
	write.table(libstats_reads,file=paste(out_prefix,".libstats_reads.tsv",sep=""),sep="\t",quote=F,row.names=F)
	write.table(libstats_pct,file=paste(out_prefix,".libstats_pct.tsv",sep=""),sep="\t",quote=F,row.names=F)
}


sc_dat$sample <- factor(sc_dat$sample, levels = unique(sc_dat$sample))
#sc_dat = read.csv(data_file,header = T,sep = "\t",stringsAsFactors = F)
#sc_dat$sample <- gsub(".*/","",sc_dat$sample)
#sc_dat$sample <- gsub(".coutc.csv$","",sc_dat$sample)
#sc_dat$sample <- factor(sc_dat$sample, levels = unique(sc_dat$sample))


if (!is.null(cell_names_path)){
  cell_names <- read.csv(cell_names_path, header=T, sep="\t", stringsAsFactors = F)
  sc_dat <- merge(sc_dat,cell_names[,c("sample","cell_idx","library","plate")],by.x = c("sample","cell_idx"),by.y = c("sample","cell_idx"))
} #else {
  #cell_names_tmp <- rbind(cell_names_tmp,data.frame(sample=dat$name[i],plate=((i-1) %/% libs_per_plate)+1, library = ((i-1) %% libs_per_plate)+1, cell_idx=seq(1+cell_idx_offset,ncol(tmp)-1+cell_idx_offset),cell_name=colnames(tmp)[2:ncol(tmp)]))
#}


str(sc_dat)

## simple way to adjust height of plots to the number of samples  
num_samples=length(levels(factor(sc_dat$sample)))
print(num_samples)
height=6;

if (num_samples>4){
	height <- height + as.integer(max(1,as.integer(num_samples/2)-2) * 1);
}

libs_per_plate = 2
if (is_split_library) libs_per_plate = 4


p<-ggplot(dat=sc_dat,aes(x=(READS_UNIQFEAT),y=(UMI),color=sample))+geom_point(size=1,alpha=0.8) + 
	facet_wrap(~sample,ncol=libs_per_plate)+
	xlab("reads on feature per cell")+
	ylab("unique UMIs per cell (trancripts)")+
	theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))

if (!is.null(plot_format)){
	ggsave(plot=p,filename=paste(out_prefix,".reads_UMI_plot.",plot_format,sep=""),device=plot_format,dpi=200,height=height)
}


## input:
## sample                cell  cell_reads  sample2               cell2  cell_transcripts
## B6_11wks_1.coutc.csv  1     40809       B6_11wks_1.coutb.csv  1      8716
## B6_11wks_1.coutc.csv  2     17736       B6_11wks_1.coutb.csv  2      4615

sc_dat <- sc_dat %>% group_by(sample) %>% mutate(sample_reads = sum(READS_UNIQFEAT)) %>% ungroup() %>% 
		mutate(cRPM=READS_UNIQFEAT/sample_reads*1000000, cUPM=UMI/sample_reads*1000000) %>% 
		mutate(cRPM_zscore = (log2(cRPM+5)-mean(log2(cRPM+5),trim=0.01))/sd(log2(cRPM+5)),
			   cUPM_zscore = (log2(cUPM+5)-mean(log2(cUPM+5),trim=0.01))/sd(log2(cUPM+5)),
			   y=-((cell_idx-1)%%16),x=((cell_idx-1)%/%16)+1)

p <- ggplot(sc_dat,aes(x=x,y=y,fill=cUPM_zscore))+ 
		geom_tile() + 
		facet_wrap(~sample,ncol = libs_per_plate, scales = "free") + 
    	#facet_grid(plate~library,scales = "free") + 
		scale_fill_gradient2(low="red",mid="blue",limits=c(-3,3),high="cyan" ) + 
		#coord_fixed() + 
    	scale_y_continuous(breaks=-seq(1,15,2),labels = as.character(seq(2,16,2))) +
    	scale_x_continuous(breaks=seq(2,(max(sc_dat$cell_idx-1)%/%16)+1,2),labels = as.character(seq(2,(max(sc_dat$cell_idx-1)%/%16)+1,2))) +
		theme_minimal() +
		theme(strip.text.x = element_text(size = 16, colour = "black",face="bold")) +
		ggtitle(paste("cell z-score of norm. transcripts per cell (cUPM)")) + 
		theme(plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5))

if (!is.null(plot_format)){
	ggsave(plot=p,filename=paste(out_prefix,".plate_cUPM.",plot_format,sep=""),device=plot_format,dpi=200,height=height)
}


p<-ggplot(sc_dat,aes(x=x,y=y,fill=cRPM_zscore))+ 
		geom_tile() + 
		facet_wrap(~sample,ncol = libs_per_plate,scales = "free") + 
		scale_fill_gradient2(low="red",mid="blue",limits=c(-3,3),high="cyan" ) + 
    	scale_y_continuous(breaks=-seq(1,15,2),labels = as.character(seq(2,16,2))) +
    	scale_x_continuous(breaks=seq(2,(max(sc_dat$cell_idx-1)%/%16)+1,2),labels = as.character(seq(2,(max(sc_dat$cell_idx-1)%/%16)+1,2))) +
		theme_minimal() +
		theme(strip.text.x = element_text(size = 16, colour = "black",face="bold")) +
		ggtitle(paste("cell z-score of norm. reads per cell (cRPM)")) + 
		theme(plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5))

if (!is.null(plot_format)){
	ggsave(plot=p,filename=paste(out_prefix,".plate_cRPM.",plot_format,sep=""),device=plot_format,dpi=200,height=height)
}

p<-ggplot(sc_dat,aes(x=x,y=y,fill=UMI))+ 
		geom_tile() + 
		facet_wrap(~sample,ncol = libs_per_plate) + 
		#scale_fill_gradient2(low="red",mid="blue",high="cyan",limits=c(min(sc_dat$cell_transcripts),max(sc_dat$cell_transcripts)),midpoint=mean(sc_dat$cell_transcripts,trim=0.05)) + 
		scale_fill_gradientn(colors=c("red","blue","cyan"),
				values=rescale(c(0,median(sc_dat$UMI)*1.00,quantile(sc_dat$UMI,probs=0.95))),
				limits=c(0,quantile(sc_dat$UMI,probs=0.95)),space = "Lab",oob=squish) +
    	scale_y_continuous(breaks=-seq(1,15,2),labels = as.character(seq(2,16,2))) +
    	scale_x_continuous(breaks=seq(2,(max(sc_dat$cell_idx-1)%/%16)+1,2),labels = as.character(seq(2,(max(sc_dat$cell_idx-1)%/%16)+1,2))) +
    	coord_fixed(ratio = 0.5) +
		theme_minimal() + 
		theme(	strip.text.x = element_text(size = 8, colour = "white",face="bold"),
				axis.text = element_text(size = 8, colour = "grey"),
				axis.ticks.y = element_blank(),
				plot.background = element_rect(fill = "black"),
				panel.grid = element_blank() ) +
        ggtitle(paste("Total number of transcripts per cell \n","red color < median (n=",median(sc_dat$UMI),") > blue clolor",sep = "")) + 
		theme(plot.title = element_text(color="white", size=10, face="bold",hjust = 0.5)) + 
		theme(legend.title = element_text(colour="grey90", size=12, face="bold")) + 
		theme(legend.text = element_text(colour="grey90", size = 12, face = "bold"))

if (!is.null(plot_format)){
	ggsave(plot=p,filename=paste(out_prefix,".plate_abs_transcripts.",plot_format,sep=""),device=plot_format,dpi=200,height=height)
}


#barcode2plate <- function(y) {
#	y1 = y-1 %% 16
#	print(y1)
#	if (y %/% 8 == 0) {
#		y2 = 2*y1-1 
#	} else {
#		y2 = ((y+1)-8) * 2
#	}
#	return(y2)
#}
