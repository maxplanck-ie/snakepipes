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


###############################################################################

files <- list.files(cellsum_path, pattern=".cellsum")

if (length(files) == 0){
	print(paste("no csv files found under",cellsum_path))
	stop(1)
}

for (f in files) {
	dat <- read.csv(f, header=T, sep="\t", stringsAsFactors = F)
	
	if (!all(c("sample","cell_idx","READS_UNIQFEAT","UMI") %in% names(dat))) {
		print("required columns not found!")
		stop(1)
	}
	
	sc_dat <- rbind(sc_dat, dat)
}
print(str(sc_dat))

sc_dat$sample <- factor(sc_dat$sample, levels = unique(sc_dat$sample))
#sc_dat = read.csv(data_file,header = T,sep = "\t",stringsAsFactors = F)

#sc_dat$sample <- gsub(".*/","",sc_dat$sample)
#sc_dat$sample <- gsub(".coutc.csv$","",sc_dat$sample)

#sc_dat$sample <- factor(sc_dat$sample, levels = unique(sc_dat$sample))

png(file=paste(out_prefix,".reads_UMI_plot.png",sep=""),width=1500,height=1500)
ggplot(dat=sc_dat,aes(x=(READS_UNIQFEAT),y=(UMI),color=sample))+geom_point(size=3,alpha=0.8) + 
	facet_wrap(~sample,ncol=4)+
	xlab("reads on feature per cell")+
	ylab("unique UMIs per cell (trancripts)")+
	theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))
dev.off()

## input:
## sample                cell  cell_reads  sample2               cell2  cell_transcripts
## B6_11wks_1.coutc.csv  1     40809       B6_11wks_1.coutb.csv  1      8716
## B6_11wks_1.coutc.csv  2     17736       B6_11wks_1.coutb.csv  2      4615

sc_dat <- sc_dat %>% group_by(sample) %>% mutate(sample_reads = sum(READS_UNIQFEAT)) %>% ungroup() %>% 
		mutate(cRPM=READS_UNIQFEAT/sample_reads*1000000, cUPM=UMI/sample_reads*1000000) %>% 
		mutate(cRPM_zscore = (log2(cRPM+5)-mean(log2(cRPM+5),trim=0.01))/sd(log2(cRPM+5)),
			   cUPM_zscore = (log2(cUPM+5)-mean(log2(cUPM+5),trim=0.01))/sd(log2(cUPM+5)),
			   y=-((cell_idx-1)%%16),x=((cell_idx-1)%/%16)+1)

num_samples=length(levels(factor(sc_dat$sample)))
print(num_samples)
height=1000;

if (num_samples>4){
	height <- height + (num_samples-4)*50;
}

print(height)

png(file=paste(out_prefix,".plate_cUPM.png",sep=""),width=1000,height=height)

ggplot(sc_dat,aes(x=x,y=y,fill=cUPM_zscore))+ 
		geom_tile() + 
		facet_wrap(~sample,ncol = 4,scales = "free") + 
		scale_fill_gradient2(low="red",mid="blue",limits=c(-3,3),high="cyan" ) + 
		coord_fixed() + 
		theme_minimal() +
		theme(strip.text.x = element_text(size = 18, colour = "black",face="bold")) +
		ggtitle(paste("cell z-score of norm. transcripts per cell (cUPM)")) + 
		theme(plot.title = element_text(color="black", size=22, face="bold",hjust = 0.5))

		
dev.off()

png(file=paste(out_prefix,".plate_cRPM.png",sep=""),width=1000,height=height)

ggplot(sc_dat,aes(x=x,y=y,fill=cRPM_zscore))+ 
		geom_tile() + 
		facet_wrap(~sample,ncol = 4,scales = "free") + 
		scale_fill_gradient2(low="red",mid="blue",limits=c(-3,3),high="cyan" ) + 
		coord_fixed() + 
		theme_minimal() +
		theme(strip.text.x = element_text(size = 18, colour = "black",face="bold")) +
		ggtitle(paste("cell z-score of norm. reads per cell (cRPM)")) + 
		theme(plot.title = element_text(color="black", size=22, face="bold",hjust = 0.5))

dev.off()


png(file=paste(out_prefix,".plate_abs_transcripts.png",sep=""),width=1000,height=height)

ggplot(sc_dat,aes(x=x,y=y,fill=cell_transcripts))+ 
		geom_tile() + 
		facet_wrap(~sample,ncol = 4,scales = "free") + 
		#scale_fill_gradient2(low="red",mid="blue",high="cyan",limits=c(min(sc_dat$cell_transcripts),max(sc_dat$cell_transcripts)),midpoint=mean(sc_dat$cell_transcripts,trim=0.05)) + 
		scale_fill_gradientn(colors=c("red","blue","cyan"),
				values=rescale(c(0,median(sc_dat$cell_transcripts),max(sc_dat$cell_transcripts))),
				limits=c(0,max(sc_dat$cell_transcripts)),space = "Lab") +
		coord_fixed() + 
		theme_minimal() + 
		theme(	strip.text.x = element_text(size = 17, colour = "white",face="bold"),
				axis.text = element_text(size = 10, colour = "grey"),
				axis.ticks.y = element_blank(),
				plot.background = element_rect(fill = "black"),
				panel.grid = element_blank() ) +
        ggtitle(paste("Total number of transcripts per cell \n","red color < median (n=",median(sc_dat$cell_transcripts),") > blue clolor",sep = "")) + 
		theme(plot.title = element_text(color="white", size=22, face="bold",hjust = 0.5)) + 
		theme(legend.title = element_text(colour="grey90", size=16, face="bold")) + 
		theme(legend.text = element_text(colour="grey90", size = 16, face = "bold"))

dev.off()

