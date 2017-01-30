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
data_file = ""
if (!is.na(args[1])) {
	data_file <- args[1] 
} else {
	print("Data path param missing! Exit...\n")
	quit(1)
}

print(data_file)

sc_dat = read.csv(data_file,header = T,sep = "\t",stringsAsFactors = F)

sc_dat$sample <- gsub(".*/","",sc_dat$sample)
sc_dat$sample <- gsub(".coutc.csv$","",sc_dat$sample)

sc_dat$sample <- factor(sc_dat$sample, levels = unique(sc_dat$sample))

png(file=paste(data_file,".reads_UMI_plot.png",sep=""),width=1500,height=1500)
ggplot(dat=sc_dat,aes(x=(cell_reads),y=(cell_transcripts),color=sample))+geom_point(size=3,alpha=0.8) + 
	facet_wrap(~sample,ncol=4)+
	xlab("reads on feature per cell")+
	ylab("unique UMIs per cell (trancripts)")+
	theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))
dev.off()

## input:
## sample                cell  cell_reads  sample2               cell2  cell_transcripts
## B6_11wks_1.coutc.csv  1     40809       B6_11wks_1.coutb.csv  1      8716
## B6_11wks_1.coutc.csv  2     17736       B6_11wks_1.coutb.csv  2      4615

sc_dat <- sc_dat %>% group_by(sample) %>% mutate(sample_reads = sum(cell_reads)) %>% ungroup() %>% 
		mutate(cRPM=cell_reads/sample_reads*1000000, cTPM=cell_transcripts/sample_reads*1000000) %>% 
		mutate(cRPM_zscore = (log2(cRPM+5)-mean(log2(cRPM+5),trim=0.01))/sd(log2(cRPM+5)),
			   cTPM_zscore = (log2(cTPM+5)-mean(log2(cTPM+5),trim=0.01))/sd(log2(cTPM+5)),
			   y=-((cell-1)%%16),x=((cell-1)%/%16)+1)

num_samples=length(levels(factor(sc_dat$sample)))
print(num_samples)
height=1000;

if (num_samples>4){
	height <- height + (num_samples-4)*50;
}

print(height)

png(file=paste(data_file,".plate_cTPM.png",sep=""),width=1000,height=height)

ggplot(sc_dat,aes(x=x,y=y,fill=cTPM_zscore))+ 
		geom_tile() + 
		facet_wrap(~sample,ncol = 4,scales = "free") + 
		scale_fill_gradient2(low="red",mid="blue",limits=c(-3,3),high="cyan" ) + 
		coord_fixed() + 
		theme_minimal() +
		theme(strip.text.x = element_text(size = 18, colour = "black",face="bold")) +
		ggtitle(paste("cell z-score of norm. transcripts per cell (cTPM)")) + 
		theme(plot.title = element_text(color="black", size=22, face="bold",hjust = 0.5))

		
dev.off()

png(file=paste(data_file,".plate_cRPM.png",sep=""),width=1000,height=height)

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


png(file=paste(data_file,".plate_abs_transcripts.png",sep=""),width=1000,height=height)

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

