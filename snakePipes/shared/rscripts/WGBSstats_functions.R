##to be called by WGBS stats Rscripts

print_sessionInfo<-function(mytext){
    message(mytext)
    sink("sessionInfo.txt")
    sessionInfo()
    sink()
    return()
}



