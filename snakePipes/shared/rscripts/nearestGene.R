#!/usr/bin/env Rscript
## adds gene ID and gene symbol annotation to bedtools closest result
.libPaths(R.home("library"))


input_bed <- snakemake@input[["bed"]]  
t2g <- snakemake@input[["t2g"]] 
gene_symbol <- snakemake@input[["gene_symbol"]]

pipeline<-snakemake@params[["pipeline"]]

output_bed<-snakemake@output[["annotated_bed"]]

flist<-list(input_bed,t2g,gene_symbol)
size_v<-unlist(lapply(flist,function(X)file.info(X)$size))

if(any(is.na(size_v),sum(size_v==0)>0)){message('Some of the input files are nonexistent or empty!')

       system(paste0('touch ',output_bed))

    }else{

    ibed_tab<-data.table::fread(input_bed,header=FALSE)
    t2g_tab<-data.table::fread(t2g,header=FALSE)
    gs_tab<-data.table::fread(gene_symbol,header=FALSE)

    if(pipeline %in% c("chipseq","ATACseq")){
        ibed_tab$GeneID<-t2g_tab$V2[match(ibed_tab$V22,t2g_tab$V1)]
        ibed_tab$GeneSymbol<-gs_tab$V2[match(ibed_tab$GeneID,gs_tab$V1)]
        obed_tab<-unique(subset(ibed_tab,select=c(paste0("V",c(1:18,23,24)),"GeneID","GeneSymbol")))
        colnames(obed_tab)<-c("Chromosome","Start","End","Width","Strand","Score","nWindows","logFC.up","logFC.down","PValue","FDR","direction","rep.test","rep.logFC","best.logFC","best.test","best.start","Name","GeneStrand","Distance","GeneID","GeneSymbol")
    }


    write.table(obed_tab,file=output_bed,row.names=FALSE,sep="\t",quote=FALSE)
}

sink(file.path(snakemake@params[["wdir"]],"nearestGene.session_info.txt"))
sessionInfo()
sink()

print("DONE..!")


