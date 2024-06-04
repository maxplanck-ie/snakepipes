#### ~~~~ Functions to Run DESeq as part of SNAKEMAKE pipeline ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Check if names of the setup table are subset of the count matrix column names
#'
#' @param countdata data.frame with featurecounts output table
#' @param _sampleSheet sample info data frame from the samplesheet
#' @param alleleSpecific TRUE/FALSE is the design allele-specific?
#' @param salmon_dir directory path for salmon output files
#' @param tx2gene_annot transcript ID to gene ID annotation
#'
#' @return NA
#' @export
#'
#' @examples
#'

.libPaths(R.home("library"))

checktable <- function(countdata = NA, sampleSheet = NA, alleleSpecific = FALSE, salmon_dir = NA, tx2gene_annot = NA) {

  ## check whether colnames are allele-specific
  print(paste0("Allele-specific counts? : ", alleleSpecific))
  if(alleleSpecific) {
    coln <- gsub("(.*)_(all|genome[1|2])", "\\1" , colnames(countdata) )
  } else {
    if(is.na(salmon_dir)) {
      coln <- colnames(countdata)
    }
  }
  ## check files
  if(!is.na(salmon_dir)) {

    #mode = Salmon : check whether salmon output files exist in Salmon dir
    if(alleleSpecific){
        files <- c(paste0(salmon_dir,"/",sampleSheet$name,".genome1.quant.sf"),paste0(salmon_dir,"/",sampleSheet$name,".genome2.quant.sf"))
    }else{
        files <- paste0(salmon_dir,"/",sampleSheet$name,".quant.sf")
    }
    #files<-dir(salmon_dir,pattern=".quant.sf",full.names=TRUE)
    filecheck <- file.exists(files)
    if(!(all(filecheck == TRUE))) {
      cat("Error! The following File(s) don't exist : ")
      cat(paste(files[filecheck == FALSE], sep = "\n") )
      quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
    } else {
      print("Importing transcript abundances")
      # import counts
      countdata <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene_annot)
    }

  } else {
    # mode = Normal : check whether sample names in the samplesheet are also present in count table header
    if ( !all( is.element(sort(sampleSheet[,1]), sort(coln)) )) {
      cat("Error! Count table column names and setup table names do NOT match!\n")
      print("Count table : ")
      print(coln)
      print("Setup table : ")
      print(as.character(sampleSheet[,1]))
      quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
    } else {
        if(alleleSpecific) {
            coln_allelic <- paste(rep(sampleSheet$name, each  = 3), c("all","genome1", "genome2"), sep = "_" )
            countdata <- countdata[,coln_allelic]

        } else {
            countdata <- countdata[,sampleSheet$name]
        }
    }
  }
  is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if( is.na(salmon_dir)){
      if(all(is.wholenumber(countdata))){
            print("All countdata is integer.")
        }else{
            print("Non-integer counts detected. The data will be rounded, as this is well within the expected sampling variation of a technical replicate.")
            countdata<-round(countdata) }
  }
  return(countdata)
}

#' DEseq basic
#'
#' @param countdata count file from featurecounts
#' @param coldata sampleSheet file
#' @param fdr fdr cutoff for DE
#'
#' @return A list of dds and ddr objects
#' @export
#'
#' @examples
#'
#'

DESeq_basic <- function(countdata, coldata, fdr, alleleSpecific = FALSE, from_salmon = FALSE, size_factors=NA, customFormula=NA) {
    cnames.sub<-unique(colnames(coldata)[2:which(colnames(coldata) %in% "condition")])
    
    if(is.na(customFormula)|customFormula==""){
      d<-as.formula(noquote(paste0("~",paste(cnames.sub,collapse="+"))))
    } else {

      d<-as.formula(paste0("~",customFormula))
    }
    

    # Normal DESeq
    print("Performing basic DESeq: test vs control")
    if(isTRUE(from_salmon)) {
      print("Using input from tximport")
        dds <- DESeq2::DESeqDataSetFromTximport(countdata,
                                  colData = coldata, design =d)
                
      } else {
          if(isTRUE(alleleSpecific)) {
            rnasamp <- dplyr::select(countdata, dplyr::ends_with("_all"))
            rownames(coldata)<-colnames(rnasamp)
            countdata<-rnasamp
           }

           dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata,
                                    colData = coldata, design =d)

      }

    if(length(size_factors) > 1) {
        print("applying size factors")
        print(size_factors)
        sizeFactors(dds) = size_factors
    }
    dds <- DESeq2::DESeq(dds)
    ddr <- DESeq2::results(dds, alpha = fdr)
    a <- DESeq2::resultsNames(dds)
    auto_coef <- a[length(a)]
    ddr_shrunk <- DESeq2::lfcShrink(dds,coef=auto_coef,type="apeglm",res=ddr)
    output <- list(dds = dds, ddr = ddr, ddr_shrunk = ddr_shrunk)
    return(output)

}

#' DESeq using allele-specific design
#'
#' @param countdata count file from featurecounts
#' @param coldata sampleSheet file
#' @param fdr fdr cutoff for DE
#'
#' @return A list of dds and ddr objects
#' @export
#'
#' @examples
#'
#'

DESeq_allelic <- function(countdata, coldata, fdr, from_salmon=FALSE, customFormula=NA) {

    # AlleleSpecific DEseq
    print("Performing Allele-specific DESeq using Interaction design : Genome2 vs Genome1")
    if(isTRUE(from_salmon)) {
        # create alleleSpecific design matrix
        coldata_allelic <- data.frame(name = colnames(as.data.frame(countdata$counts)),
                   allele = rep(c("genome1", "genome2"), nrow(coldata)),
                   condition = rep(coldata$condition, each = 2) )
        rownames(coldata_allelic)<-colnames(as.data.frame(countdata$counts))
        coldata_allelic$allele<-factor(coldata_allelic$allele,levels=c("genome1","genome2"))
        coldata_allelic$condition<-factor(coldata_allelic$condition,levels=unique(coldata_allelic$condition))
        print("Using input from tximport")
        dds <- DESeq2::DESeqDataSetFromTximport(countdata,
                                  colData = coldata_allelic, design =~1)

    } else {
      print("Using input from count table")
    rnasamp <- dplyr::select(countdata, -dplyr::ends_with("_all"))

    # create alleleSpecific design matrix
    coldata_allelic <- data.frame(name = colnames(rnasamp),
                   allele = rep(c("genome1", "genome2"), nrow(coldata)),
                   condition = rep(coldata$condition, each = 2) )
    rownames(coldata_allelic)<-colnames(rnasamp)
    coldata_allelic$allele<-factor(coldata_allelic$allele,levels=c("genome1","genome2"))
    coldata_allelic$condition<-factor(coldata_allelic$condition,levels=unique(coldata_allelic$condition))
    dds <- DESeq2::DESeqDataSetFromMatrix(rnasamp, colData = coldata_allelic,
                              design = ~1)
    rownames(dds) <- rownames(rnasamp)

    }
    
    # Run DESeq
    if(length(unique(coldata_allelic$condition))>1){
      DESeq2::design(dds) <- formula(~allele + condition + allele:condition)
      dds <- DESeq2::DESeq(dds,betaPrior = FALSE)
      ddr <- DESeq2::results(dds, name=paste0("allelegenome2.condition",unique(coldata$condition)[2]))
      ddr_shrunk <- DESeq2::lfcShrink(dds,coef=paste0("allelegenome2.condition",unique(coldata$condition)[2]),type="apeglm",res=ddr)
    } else {
        DESeq2::design(dds) <- formula(~allele)
        dds <- DESeq2::DESeq(dds,betaPrior = FALSE)
        ddr <- DESeq2::results(dds, name="allele_genome2_vs_genome1")
        ddr_shrunk <- DESeq2::lfcShrink(dds,coef="allele_genome2_vs_genome1",type="apeglm",res=ddr)
    }
    output <- list(dds = dds, ddr = ddr, ddr_shrunk=ddr_shrunk)
    return(output)
}


#' Plotting and annotation of DESeq outputs
#'
#' @param DEseqout output from DEseq_basic or DEseq_allelic wrapper
#' @param countdata count file from featurecounts
#' @param coldata sampleSheet file
#' @param fdr fdr cutoff for DE
#' @param outprefix prefix of output plots and files
#' @param heatmap_topN top N genes to plot in the heatmap
#' @param geneNamesFile tsv file with gene symbols corresponding to ensembl gene IDs
#'
#' @return DEseq results (.tsv) and plots (.pdf)
#' @export
#' @import ggplot2
#'
#' @examples
#'
#'

DESeq_writeOutput <- function(DEseqout,
                              fdr,
                              geneNamesFile,
                              outprefix) {

    # get dds and ddr
    dds <- DEseqout$dds
    ddr <- DEseqout$ddr
    ddr_shrunk <- DEseqout$ddr_shrunk
    # Now prep files for outputs
    ddr.df <- as.data.frame(ddr)
    print("Adding Status: UP/DOWN")
    # Add UP/DOWN status to ddr.df
    ddr.df$Status <- ifelse(is.na(ddr.df$padj), "None",
                    ifelse(ddr.df$padj < fdr,
                         ifelse(ddr.df$log2FoldChange > 0, "UP", "DOWN"),
                         "None"))
    if (!is.null(ddr_shrunk)) {
        ddr_shrunk.df <- as.data.frame(ddr_shrunk)
        ddr_shrunk.df$Status <- ifelse(is.na(ddr_shrunk.df$padj), "None",
                    ifelse(ddr_shrunk.df$padj < fdr,
                         ifelse(ddr_shrunk.df$log2FoldChange > 0, "UP", "DOWN"),
                         "None"))
    }

    # If gene names given, then add this info to ddr.df
    if (file.exists(geneNamesFile)) {
        print("Gene symbols file found. Adding gene symbols")
        geneNames <- read.delim(geneNamesFile, header = F, stringsAsFactors=FALSE)
        geneNames <- geneNames[!duplicated(geneNames[,1]),]
        colnames(geneNames) <- c("GeneID", "external_gene_name")
        ddr.df <- merge(ddr.df, geneNames, by.x = 0, by.y = "GeneID" , all.x = TRUE)
        rownames(ddr.df) <- ddr.df[,1]
        ddr.df[,1] <- NULL
        if (!is.null(ddr_shrunk)) {
            ddr_shrunk.df <- merge(ddr_shrunk.df, geneNames, by.x = 0, by.y = "GeneID" , all.x = TRUE)
            rownames(ddr_shrunk.df) <- ddr_shrunk.df[,1]
            ddr_shrunk.df[,1] <- NULL
            }
    } else {
        print("Gene symbols file not found. Gene IDs used")
    }
    print("writing results")
    ## Write back diffExp output
    write.table(ddr.df,file = paste0(outprefix, "_DEresults.tsv"),sep = "\t",quote = FALSE)
    if (!is.null(ddr_shrunk)){
        write.table(ddr_shrunk.df,file = paste0(outprefix, "_DEresults_LFCshrunk.tsv"),sep = "\t",quote = FALSE)
        }
    # write normalized counts
    write.table(DESeq2::counts(dds, normalized = T),
            file = paste0(outprefix, "_counts_DESeq2.normalized.tsv"),
                      sep = "\t", quote = FALSE, col.names = NA)
    save(dds, ddr, ddr_shrunk, file = paste0(outprefix,"_DESeq.Rdata"))

}

## modified DeSeq downstream function
DESeq_downstream <- function(DEseqoutRdata,## Rdata object written by the DESeq_writeOutput function
                     countdata,
                     coldata,
                     fdr,
                     outprefix,
                     heatmap_topN,
                     geneNamesFile) {
    # get dds and ddr
    load(DEseqoutRdata)

    df.filt <- ddr.df[which(ddr.df$padj < fdr),]
    df.plot <- data.frame(Status = c("Up","Down"),
                    Genes = c(sum(df.filt$Status == "UP"),
                            sum(df.filt$Status == "DOWN"))
    )
    print("Filtered data set")
    print(head(df.filt))

    ## Rlog transform
    rld <- DESeq2::rlog(dds)

    ## Expression density data (add mean and independent filtering threshold)
    print("Preparing data: expression density")
    toplot <- data.frame(DESeq2::counts(dds, normalized = T))
    toplot <- stack(toplot, select = colnames(toplot))
    ind_filt_thres <- as.numeric(S4Vectors::metadata(ddr)$filterThreshold)
    # plotdata
    pld <- ggplot(toplot, aes(values, colour = ind, alpha = 0.5)) +
        geom_line(aes(color = ind), stat = "density", alpha = 0.5) +
        scale_x_log10(name = "\nnormalized counts",
                  breaks = c(0.1,1,10,100,1000,10000,100000),
                  limits = c(0.1,100000) ) +
        scale_y_continuous(name = "density\n") +
        scale_colour_discrete(name = "Samples") +
        geom_vline(xintercept = 10, colour = "grey", linetype = "dashed") +
        theme_minimal() +
        ggtitle("Density plot\n") +
        theme()

    ## Sample distances
    print("Preparing data: sample distances")
    sampleDists <- dist(t(SummarizedExperiment::assay(rld)))

    ## Euclidean sample distance heatmap
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- sprintf("%s\n(%s)", colnames(rld), rld$condition)
    colnames(sampleDistMatrix) <- sprintf("%s\n(%s)", colnames(rld), rld$condition)
    colours <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "GnBu")))(255)

    ## PCA data
    print("Preparing data: PCA")
    PCAdata <- DESeq2::plotPCA(rld, intgroup = c("name", "condition"), returnData=TRUE)
    percentVar <- round(100 * attr(PCAdata, "percentVar"))


    ## Take out cooks distance statistics for outlier detection
    W <- ddr$stat
    maxCooks <- apply(SummarizedExperiment::assays(dds)[["cooks"]],1,max)
    idx <- !is.na(W)
    m <- ncol(dds)
    p <- 3

    ## Extract data for heatmap
    # order by fold change (by abs foldch if only few top genes requested)
    print("Preparing data: heatmap")
    if (nrow(df.filt) > 0) {
        d <- data.frame(id = rownames(df.filt), padj = df.filt$padj)
        if (length(rownames(d)) < heatmap_topN ) {
            heatmap_topN <- nrow(d)
        }

        d_topx_padj <- d[order(d$padj, decreasing = F),][1:heatmap_topN,]
        heatmap_data <- SummarizedExperiment::assay(rld)[as.character(d_topx_padj$id),]

        ## create another df to get the gene names
        if (file.exists(geneNamesFile)) {
            mdf <- data.frame(row.names = rownames(ddr.df), external_gene_name = ddr.df$external_gene_name)
            htdat <- merge(mdf, heatmap_data, by = 0, all.y = TRUE, order = FALSE)
            napos <- which(is.na(htdat$external_gene_name))
            htdat[napos, "external_gene_name"] <- as.character(htdat[napos, "Row.names"])

            # replace rownames of heatmap data with these gene names
            rownames(heatmap_data) <- htdat$external_gene_name[match(rownames(heatmap_data),htdat$Row.names)]
        }
        ## scaling is not so useful for already normalized data, or?
        #heatmap_data <- scale(heatmap_data, center = TRUE, scale = TRUE)
    } else {
        print("No DE genes detected!")
    }

    ## Write output
    print("Plotting...")

    pdf(paste0(outprefix, "_plots.pdf"))
    # 1. sparsity
    DESeq2::plotSparsity(dds)

    # 2. expression density plot (add mean and independent filtering threshold)
    print(pld + geom_line(data = ddr.df, aes(baseMean),
                stat = "density", alpha = 0.5, colour = "blue", size = 1.1) +
        geom_vline(xintercept = ind_filt_thres, colour = "red", size=1.1) +
        ggtitle(sprintf("Density plot\n(independent filtering: %.3f)\n", ind_filt_thres)) +
        theme(legend.position = "bottom")
    )

    # 3. Dispersion estimates
    DESeq2::plotDispEsts(dds)

    # 4. Sample distances heatmap
    pheatmap::pheatmap(sampleDistMatrix, color = colours, display_numbers = TRUE)

    # 5. PCA plots
    print(ggplot(PCAdata, aes(PC1, PC2, color = condition, shape = name)) +
        geom_hline(aes(yintercept = 0), colour = "grey") +
        geom_vline(aes(xintercept = 0), colour = "grey") +
        geom_point(size = 5) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_bw(base_size = 14) +
        ggtitle("PCA\n") +
        scale_shape_manual(values = c(0:18,33:17))
    )

    # 6. Cooks distance
    plot(rank(W[idx]), maxCooks[idx], xlab = "rank of Wald statistic",
         ylab = "maximum Cook's distance per gene",
         ylim = c(0,5), cex = .4, col = rgb(0,0,0,.3))
    abline(h = qf(.99, p, m - p))

    # 7. Heatmap topN genes rld
    pheatmap::pheatmap(heatmap_data,
                 cluster_rows = TRUE,
                 clustering_method = "average",
                 show_rownames = TRUE,
                 cluster_cols = FALSE,
                 color = colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(255),
                 main = sprintf("Heatmap : Top %d DE genes (by p-value) color: norm. expression (rld) ", heatmap_topN)
    )

    # 7.1 Heatmap topN genes z-score
    pheatmap::pheatmap(heatmap_data,
            cluster_rows = TRUE,
            clustering_method = "average",
            show_rownames = TRUE,
            cluster_cols = FALSE,
            scale = "row",
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdBu")))(255),
            main = sprintf("Heatmap : Top %d DE genes (by p-value) color: z-score  ", heatmap_topN)
    )

    # 8. MAplot
    DESeq2::plotMA(ddr)

    # 9. Volcano plot
    volcano_xlim <- max(abs(quantile(na.omit(ddr.df$log2FoldChange),p= c(0.01, 0.99)) ) )
    volcano_ylim <- quantile(-log(na.omit(ddr.df$padj) ), p = 0.99)
    print(ggplot(ddr.df, aes(log2FoldChange, -log10(padj))) +
        geom_point(alpha = 0.4, size = 1, colour = "grey50") +
        scale_x_continuous(limits = c(-volcano_xlim, volcano_xlim)) +
        scale_y_continuous(limits = c(0, volcano_ylim) ) +
        geom_vline(xintercept = 0, col = "navy") +
        geom_point(data = df.filt, aes(log2FoldChange, -log10(padj)), alpha = 0.6, size = 1, colour = "red" ) +
        theme_bw(base_size = 14) +
        labs(x = "log2-fold change",
             y = "-log10 (padj)",
             title = sprintf("Volcano plot\n(FDR: %.2f, up: %d, down: %d)",
                         fdr, df.plot$Genes[df.plot$Status == 'Up'],
                         df.plot$Genes[df.plot$Status == 'Down']))
    )

    # 9. P-value distribution
    print(ggplot(ddr.df) +
        geom_histogram(aes(pvalue, fill = "a"), colour = "grey20", alpha = 0.5, stat = "bin") +
        geom_histogram(aes(padj, fill = "b"), colour = "grey20", alpha = 0.5, stat = "bin") +
        scale_fill_manual(name = "group", values = c("a" = "steelblue", "b" = "grey20"),
                    labels = c("a" = "p-value", "b" = "padj")) +
        geom_vline(xintercept = fdr, colour = "red")
    )


    # 10. UP/DOWN barplot
    print(ggplot(df.plot,ggplot2::aes(Status, Genes, fill = Status)) +
             geom_bar(stat = "identity", position = "dodge")
    )

    dev.off()
}
