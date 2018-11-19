##to be called by WGBS stats Rscripts

print_sessionInfo<-function(mytext){
    message(mytext)
    sink("sessionInfo.txt")
    print(sessionInfo())
    sink()
    }

get_nrow_input<-function(resdir,statscat){
    if(statscat %in% "single_CpGs"){
        load(file.path(resdir,"limdat.LG.RData"))
        nri<-nrow(limdat.LG)
    } else if (statscat %in% "metilene_DMRs"){
        nri<-length(readLines(file.path(resdir,"singleCpG.metilene.bed")))
    } else if (statscat %in% "user_intervals"){
        nri<-unlist(lapply(dir(resdir,pattern="*.aggCpG.RData",full.names=TRUE),function(X){load(X);nrow(bedtab)}))
    }
    return(nri)
}

get_nrow_filtered<-function(resdir,statscat){
    if(statscat %in% "single_CpGs"){
        load(file.path(resdir,"limdat.LG.RData"))
        nrf<-sum(complete.cases(limdat.LG))
    } else if (statscat %in% "metilene_DMRs"){
        load(file.path(resdir,"CGI.limdat.CC.RData"))
        nrf<-nrow(CGI.limdat.CC)
    } else if (statscat %in% "user_intervals"){
        nrf<-unlist(lapply(dir(resdir,pattern="*.aggCpG.RData",full.names=TRUE),function(X){load(X);nrow(CGI.limdat.CC)}))
    }
    return(nrf)
}

get_nrow_topTable<-function(resdir,statscat){
    if(statscat %in% "single_CpGs"){
        load(file.path(resdir,"singleCpG.RData"))
        nrtt<-nrow(tT_filt)
    } else if (statscat %in% "metilene_DMRs"){
        nrtt<-length(readLines(file.path(resdir,"singleCpG.metilene.limma_filtered.bed")))
    } else if (statscat %in% "user_intervals"){
        nrtt<-unlist(lapply(dir(resdir,pattern="*.tT_filt.txt",full.names=TRUE),function(X)length(readLines(X))))
    }
    return(nrtt)
}



