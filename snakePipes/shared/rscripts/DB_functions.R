#### ~~~~ Functions to Run CSAW as part of SNAKEMAKE pipeline ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Read the Files and Count windows for ChIPseq Samples
#'
#' @param sampleSheet tsvfile with sample information
#' @param fragmentLength fragment length of sequencing
#' @param window_size window size to count the reads in
#' @param alleleSpecific TRUE/FALSE whether samples have to be compared with allele-specific design
#' @param pe.param parameters to read bam files
#' @return chipCountObject : a list with window counts and sampleSheet
#' @export
#' @examples
#' readfiles_chip(csvFile = "testBAMs/testSampleSheet.csv", refAllele = "pat")
#'

.libPaths(R.home("library"))

readfiles_chip <- function(sampleSheet, fragmentLength, window_size, alleleSpecific = FALSE, pe.param){

    # check that not >2 conditions are given
    if(length(unique(sampleSheet$condition)) > 2 ) {
        stop("only up to 2 conditions can be used for Differential binding analysis")
    }
    if( isTRUE(alleleSpecific) ) {
        message("Mode : AlleleSpecific")
        # make allele-specific model matrix for edgeR
        design <- data.frame(name = rep(sampleSheet$name, each = 2),
                       condition = factor(rep(sampleSheet$condition, each = 2)),
                       allele = factor(rep(c("genome1", "genome2"), nrow(sampleSheet)) )
                       )
        design$allele <- relevel(design$allele, "genome1")
        design$condition <- droplevels(design$condition)

        if(length(unique(sampleSheet$condition)) == 1 ) {
            # for 1 samples, use normal design
            message("1 sample used : comparing genome2 to genome1")
            designm <- model.matrix(~ allele, data = design)
            rownames(designm)<-paste0(design$name,".",design$allele)
            designType <- "allele"
        } else {
            # for 1 sample, use interaction design
            message(">1 samples used : comparing genome2 to genome1 blocking for different conditions")
            designm <- model.matrix(~ allele + condition,data = design)
            rownames(designm)<-paste0(design$name,".",design$allele)
            designType <- "blocking"
        }

        # define bam files to read
        bam.files <- list.files("../allelic_bams",
                        pattern = paste0(sampleSheet$name,".genome[1-2].sorted.bam$", collapse = "|"),
                        full.names = TRUE )
    } else {
        message("Mode : Differential Binding")
        cnames.sub<-unique(colnames(sampleSheet)[2:which(colnames(sampleSheet) %in% "condition")])
        d<-as.formula(noquote(paste0("~",paste(cnames.sub,collapse="+"))))
        # make model matrix for differential binding
        sampleSheet$condition = factor(sampleSheet$condition )
        sampleSheet$condition <- relevel(sampleSheet$condition, ref = as.character(sampleSheet$condition[1]))# make the first entry the base level
        designm <- model.matrix(d, data = sampleSheet)
        rownames(designm)<-sampleSheet$name
        designType <- "condition"
        # define bam files to read
        bam.files <- list.files(paste0("../",bam_folder),
                                pattern = paste0(sampleSheet$name,bam_pfx,".bam$", collapse = "|"),
                                full.names = TRUE )
        
    }

    message("bam files used: ")
    message(bam.files)
    # readFiles using CSAW
    mincount <- 20
    message(paste0("Counting reads in windows.. windows with total counts < ", mincount, " are discarded"))
    counts <- csaw::windowCounts(bam.files = bam.files, param = pe.param, ext = fragmentLength, spacing = window_size, filter = mincount)
    if (designType==("condition")){
        colnames(counts)<-gsub(paste0(bam_pfx,".bam"),"",basename(bam.files))
        counts<-counts[,sampleSheet$name]} else {
        colnames(counts)<-gsub(".sorted.bam","",basename(bam.files))
        counts<-counts[,paste0(rep(sampleSheet$name,each=2),".genome",c(1,2))]
    }
    print(head(counts))

    # output
    chipCountObject <- list(windowCounts = counts, sampleSheet = sampleSheet,
                            design = designm, designType = designType, pe.param = pe.param)
    return(chipCountObject)
}

#' Make plots to select window size and pe-distance cutoffs
#'
#' @param bam.file bam file to use
#' @param outplot path to output pdf
#' @param pe.param parameters to read bam files
#' @return QC plots as pdf
#' @export
#' @examples
#' makeQCplots_chip(bam.file, outplot)
#'

makeQCplots_chip_PE <- function(bam.file, outplot, pe.param){

    ## Histogram to check frag size cutoff
    message("Checking fragment sizes")
    fragsize <- csaw::getPESizes(bam.file)

    ## Checking cross-correlation
    message("Checking strand cross-correlation")
    max.delay <- 500
    dedup.on <- csaw::readParam(dedup = TRUE, minq = 20)
    CCF <- csaw::correlateReads(bam.file, max.delay, param = dedup.on)

    ## Choosing appropriate window size
    message("Checking read distribution around putative peaks")
    plotwc <- function(curbam){
        windowed <- csaw::windowCounts(curbam, spacing = 50, param = pe.param, filter = 20)
        rwsms <- rowSums(SummarizedExperiment::assay(windowed))
        maxed <- csaw::findMaxima(SummarizedExperiment::rowRanges(windowed), range = 1000, metric = rwsms)
        curbam.out <- csaw::profileSites(curbam, SummarizedExperiment::rowRanges(windowed)[maxed], param = pe.param)
        return(curbam.out)
    }
    collected <- plotwc(bam.file)
    xranged <- as.integer(names(collected))

    ## plot
    message("Plotting")
    pdf(outplot)
    # fragment sizes
    hist(log10(fragsize$sizes),
         breaks = 50,
         xlab= " log10(Fragment sizes)",
         ylab = "Frequency",
         main = "fragment sizes",
         col = "steelblue")
    abline(v = 400, col = "red",lwd = 3)

    # cross correlation
    plot(0:max.delay, CCF, type = "l", ylab = "CCF", xlab = "Delay (bp)", main = "PE-Cross-correlation")

    # coverage in windows
    plot(xranged, collected, type = "l", col = "blue", xlim = c(-1000, 1000), lwd = 2,
         xlab = "Distance (bp)", ylab = "Relative coverage per base")
    abline(v = c(-150,200), col = "dodgerblue", lty = 2)
    legend("topright", col = "dodgerblue", legend = "specified window size")

    dev.off()

}

makeQCplots_chip_SE <- function(bam.file, outplot, pe.param){

    ## Checking cross-correlation
    message("Checking strand cross-correlation")
    max.delay <- 500
    dedup.on <- csaw::readParam(dedup = TRUE, minq = 20)
    CCF <- csaw::correlateReads(bam.file, max.delay, param = dedup.on)

    ## Choosing appropriate window size
    message("Checking read distribution around putative peaks")
    plotwc <- function(curbam){
        windowed <- csaw::windowCounts(curbam, spacing = 50, param = pe.param, filter = 20)
        rwsms <- rowSums(SummarizedExperiment::assay(windowed))
        maxed <- csaw::findMaxima(SummarizedExperiment::rowRanges(windowed), range = 1000, metric = rwsms)
        curbam.out <- csaw::profileSites(curbam, SummarizedExperiment::rowRanges(windowed)[maxed], param = pe.param)
        return(curbam.out)
    }
    collected <- plotwc(bam.file)
    xranged <- as.integer(names(collected))

    ## plot
    message("Plotting")
    pdf(outplot)
    
    # cross correlation
    plot(0:max.delay, CCF, type = "l", ylab = "CCF", xlab = "Delay (bp)", main = "SE-Cross-correlation")

    # coverage in windows
    plot(xranged, collected, type = "l", col = "blue", xlim = c(-1000, 1000), lwd = 2,
         xlab = "Distance (bp)", ylab = "Relative coverage per base")
    abline(v = c(-150,200), col = "dodgerblue", lty = 2)
    legend("topright", col = "dodgerblue", legend = "specified window size")

    dev.off()

}


#' TMM normalize (get the normfactors out) using given window size
#'
#' @param chipCountObject output from readfiles_chip command
#' @param binsize Size of bins to calculate the normalization factors
#' @param plotfile file with output plots
#' @return Normalized chipCountObject
#' @export
#' @examples
#' tmmNormalize_chip(chipCountObject,binsize = 10000, plotfile = "TMM_normalizedCounts.pdf")
#'

tmmNormalize_chip <- function(chipCountObject, binsize, plotfile){


    bam.files <- SummarizedExperiment::colData(chipCountObject$windowCounts)$bam.files
    # Get norm factors
    wider <- csaw::windowCounts(bam.files, bin = TRUE, width = binsize, param = chipCountObject$pe.param)
    if(useSpikeInForNorm){
        tab<-read.table(scale_factors,sep="\t",header=TRUE,as.is=TRUE,quote="")
        normfacs<-1/(tab$scalingFactor[match(colnames(chipCountObject$windowCounts),tab$sample)]) }else{
        normfacs <- csaw::normFactors(wider, se.out=FALSE)}

    chipCountObject$normFactors <- normfacs

    # get norm counts
    adj.counts <- edgeR::cpm(csaw::asDGEList(wider), log = TRUE)
    chipCountObject$background_logcpm <- adj.counts

    # plot normalized counts
    pdf(plotfile)
    #par(mfrow = c(3, 3), mar = c(5, 4, 2, 1.5))
    #for (i in 1:(length(bam.files) - 1)) {
    #    cur.x <- adj.counts[, 1]
    #    cur.y <- adj.counts[, 1 + i]
    #    smoothScatter(x = (cur.x + cur.y)/2 + 6*log2(10),
    #                y = cur.x-cur.y, xlab = "A",
    #                ylab = "M",
    #                main = paste("1 vs", i + 1))

    #    all.dist <- diff(log2(normfacs[c(i + 1, 1)]))
    #    abline(h = all.dist, col = "red")
    #}
    ## MDS plot to check for replicate variability
    for (top in c(100, 500, 1000, 5000)) {
        limma::plotMDS(adj.counts, main = top,
                   col = as.numeric(chipCountObject$sampleSheet$condition),
                   labels = chipCountObject$sampleSheet$name, top = top)
    }
    dev.off()

    ## Return normfactors
    return(chipCountObject)
}


#' Test for Diff Bound windows using EdgeR (then merge windows into regions)
#'
#' @param chipCountObject output from tmmNormalize_chip
#' @param plotfile file with output plots
#' @param tfname which TF to extract results for (must match with the name in sampleSheet)
#' @return chipResultObject with differentially bound regions
#' @export
#' @examples
#' getDBregions_chip(chipCountObject,plotfile = NULL, tfname = "msl2")
#'

getDBregions_chip <- function(chipCountObject, plotfile = NULL){

    # Make DGElist
    y <- csaw::asDGEList(chipCountObject$windowCounts, norm.factors = chipCountObject$normFactors)
    if(chipCountObject$designType=="condition"){
    colnames(y)<-chipCountObject$sampleSheet$name}else{
    colnames(y)<-paste0(rep(chipCountObject$sampleSheet$name,each=2),".genome",c(1,2))}
    design <- chipCountObject$design
    # Estimate dispersions
    y <- edgeR::estimateDisp(y, design)
    o <- order(y$AveLogCPM)
    fit <- edgeR::glmQLFit(y, design, robust=TRUE)

    # and plot dispersions
    if(!(is.null(plotfile))){
        pdf(plotfile)
        par(mfrow = c(1,2))
        plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type = "l", lwd = 2,
             ylim = c(0, 1), xlab = expression("Ave."~Log[2]~"CPM"),
             ylab = ("Biological coefficient of variation"))
        edgeR::plotQLDisp(fit)
        dev.off()
    }

    ### TEST for DB windows
    # check design type
    if(chipCountObject$designType != "condition") {
        results <- edgeR::glmQLFTest(fit, coef = paste0("allelegenome2"))
    } else {
        results <- edgeR::glmQLFTest(fit, coef = paste0("condition",unique(chipCountObject$sampleSheet$condition)[2]))
    }

    # Merge DB windows into regions: Using quick and dirty method
    merged <- csaw::mergeWindows(SummarizedExperiment::rowRanges(chipCountObject$windowCounts), tol = 100L)
    # get combined test p-value for merged windows
    tabcom <- csaw::combineTests(merged$id, results$table, pval.col = 4, fc.col = 1)
    # get fold change of the best window within each combined cluster
    tab.best <- csaw::getBestTest(merged$id, results$table,pval.col=4,cpm.col=1)
    tabcom$best.logFC <- tab.best$rep.logFC
    tabcom$best.test <- tab.best$rep.test
    tabcom$best.start <- GenomicRanges::start(SummarizedExperiment::rowRanges(chipCountObject$windowCounts))[tab.best$rep.test]

    # Return all results
    chipResultObject <- list(fit = fit, results = results, mergedRegions = merged, combinedPvalues = tabcom)
    return(chipResultObject)
}


#' Annotate and print the output regions
#'
#' @param chipResultObject output from getDBregions_chip
#' @param outfile_prefix name of output files
#' @param fdrcutoff fdr cutoff to call significantly bound regions
#' @return File with differentially bound regions
#'
#' @export
#' @examples
#' writeOutput_chip(chipResultObject, outfile_prefix)
#'

writeOutput_chip <- function(chipResultObject, outfile_prefix, fdrcutoff,lfccutoff){
    # get merged regions
    merged <- chipResultObject$mergedRegions
    tabcom <- chipResultObject$combinedPvalues
    merged$region$score <- -10*log10(tabcom$FDR)
    names(merged$region) <- paste0("region", 1:length(merged$region))
    tabcom$name <- names(merged$region)
    ## export merged data
    rtracklayer::export.bed(merged$region, paste0(outfile_prefix, "_allregions.bed"))
    write.table(tabcom, file = paste0(outfile_prefix,"_scores.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

    ## export FDR significant data
    is.sig <- tabcom$FDR <= fdrcutoff
    test <- merged$region[is.sig]

    if(length(test) > 0){
        rtracklayer::export.bed(test, paste0(outfile_prefix,"_significant.bed"))
    } else {
        warning("output empty! please lower the fdr threshold.")
    }
    ##merge regions with stats
    print(head(as.data.frame(merged$region)))
    print(head(tabcom))
    tabx<-as.data.frame(merged$region,stringsAsFactors=FALSE)
    tabx$name<-rownames(tabx)
    full_res<-as.data.frame(merge(x=tabx,y=tabcom,by.x="name",by.y="name"),stringsAsFactors=FALSE) 
    full_res<-full_res[,c(2:ncol(full_res),1)]
    print(sprintf("Colnames of result file are %s",colnames(full_res)))
    full_res[,2]<-full_res[,2]-1
    full_res[,2]<-format(full_res[,2], scientific = FALSE,trim=TRUE)
    full_res[,3]<-format(full_res[,3], scientific = FALSE,trim=TRUE)
    ##filter full result for FDR and LFC and write to output
    full_res.filt<-subset(full_res,(FDR<=fdrcutoff)&(abs(best.logFC)>=lfccutoff))
    if(nrow(full_res.filt)>0){
    write.table(full_res.filt,file="Filtered.results.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.bed")}
    res.filt.up<-subset(full_res.filt,direction %in% "up")
    if(nrow(res.filt.up)>0){
    write.table(res.filt.up,file="Filtered.results.UP.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.UP.bed")}
    res.filt.down<-subset(full_res.filt,direction %in% "down")
    if(nrow(res.filt.down)>0){
    write.table(res.filt.down,file="Filtered.results.DOWN.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.DOWN.bed")}
    res.filt.mixed<-subset(full_res.filt,direction %in% "mixed")
    if(nrow(res.filt.mixed)>0){
    write.table(res.filt.mixed,file="Filtered.results.MIXED.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.MIXED.bed")}
}




######################################      UNUSED FUNCTIONS    ##########################################

#' Filtering using global bg
#'
#' @param chipCountObject output from readfiles_chip
#' @return Filtered chipCountObject
#' @export
#' @examples
#' filterByInput_chip(chipCountObject,priorCount = 5)
#'

filterByGlobal_chip <- function(chipCountObject){
    bin.size <- 2000L
    countdat <- chipCountObject$windowCounts
    binned <- windowCounts(SummarizedExperiment::colData(countdat)$bam.files, bin = TRUE,
                     width = bin.size, param = chipCountObject$pe.param)

    filter.stat <- filterWindows(countdat, background = binned, type = "global")
    keep <- filter.stat$filter > log2(3)
    print(sum(keep))
    chipCountObject$windowCounts <-  countdat[keep,]
    return(chipCountObject)
}


#' Filtering using Input windows
#'
#' @param chipCountObject output from readfiles_chip
#' @param controlbam control bam file
#' @param chipbam chip bam file
#' @param priorCount Minimum count cutoff for windows
#' @return Filtered chipCountObject
#' @export
#' @examples
#' filterByInput_chip(chipCountObject,priorCount = 5)
#'

filterByInput_chip <- function(chipCountObject, controlbam, chipbam, priorCount = 5){
    countdat <- chipCountObject$windowCounts
    control <- countdat[,which(SummarizedExperiment::colData(countdat)$bam.files %in% controlbam)]
    chip <- countdat[,which(SummarizedExperiment::colData(countdat)$bam.files %in% chipbam)]
    # Filter chip by control counts
    filter.stat <- csaw::filterWindows(chip, control, type = "control", prior.count = priorCount) # min count in window should be 5
    keep <- filter.stat$filter > log2(3)
    print(sum(keep))
    countdat <- countdat[keep,] # now "countdat" contains both input and chip counts, but filtered

    # return the same chipCountObject back, but this time filtered
    chipCountObject$windowCounts <- countdat
    return(chipCountObject)
}
