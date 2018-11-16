##to be called by WGBS stats Rscripts

print_sessionInfo<-function(mytext){
    message(mytext)
    sink("sessionInfo.txt")
    print(sessionInfo())
    sink()
    }

get_nrow_topTable<-function(resdir,statscat){
    if(statscat %in% "single_CpGs"){
        load(file.path(resdir,"singleCpG.RData"))
        nrtt<-nrow(tT_filt)
    } else if (statscat %in% "metilene_DMRs"){
        nrtt<-system(paste0('wc -l ',file.path(resdir,"singleCpG.metilene.limma_filtered.bed")),intern=TRUE)
    } else if (statscat %in% "user_intervals"){
        nrtt<-unlist(lapply(dir(resdir,pattern="*.tT_filt.txt",full.names=TRUE),function(X)system(paste0('wc -l ',X),intern=TRUE)))
    }
    return(nrtt)
}



