.libPaths(R.home("library"))

#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

importfunc<-commandArgs(trailingOnly=TRUE)[4]
source(importfunc)

options(stringsAsFactors=FALSE,na.rm=TRUE)

spath<-commandArgs(trailingOnly=TRUE)[2]
sampleSheet<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)
if(!"PlottingID" %in% colnames(sampleSheet)){sampleSheet$PlottingID<-sampleSheet$name}

############read in methylation tracks; use (reference)-sorted bed files #############################################
#sorted bed file was used at the beginning to extract CGs from intervals of interest + 0-padding was used -> all files should have same row ordering (although possibly different from the reference)
###deal with non-0-padded tables from methylDackel and BisSNP

mpath<-commandArgs(trailingOnly=TRUE)[3]
mdir<-dir(mpath,pattern="*CpG.filt2.bed",full.names=TRUE)
mshort<-gsub(".CpG.filt2.bed","",basename(mdir))
mdir<-mdir[match(sampleSheet$name,mshort)]
mshort<-gsub(".CpG.filt2.bed","",basename(mdir))


require(data.table)

mlist<-vector("list",length(mdir))
for(i in seq_along(mdir)){
    tabi<-fread(mdir[i],select=c("Beta","ms"),sep="\t",header=TRUE)
    colnames(tabi)[colnames(tabi) %in% "Beta"]<-mshort[i]
    mlist[[i]]<-tabi
}

limdat<-Reduce(function(...) merge(..., all=T,by="ms",sort=FALSE), mlist)
limdat<-limdat[,c(1,match(sampleSheet$name,colnames(limdat))),with=FALSE]


limdat.LG<-limdat
limdat.LG[,2:ncol(limdat.LG)]<-limdat.LG[,2:ncol(limdat.LG)]/100

save(limdat.LG,file="limdat.LG.RData")

limdat.LG.CC<-limdat.LG[complete.cases(limdat.LG),] 
if(nrow(limdat.LG.CC)==0){ print_sessionInfo("None of the single CpG sites passed the filtering.")}else{

    #### prepare metilene input

    limdat.LG.CC.tw<-limdat.LG.CC

    ##reorder input data so that Treatment or WT go first
    if ("Mut" %in% sampleSheet$condition){
        limdat.LG.CC.tw<-limdat.LG.CC[,c("ms",colnames(limdat.LG.CC)[match(sampleSheet$name[sampleSheet$condition %in% "Mut"],colnames(limdat.LG.CC))],colnames(limdat.LG.CC)[match(sampleSheet$name[!sampleSheet$condition %in% "Mut"],colnames(limdat.LG.CC))]),with=FALSE]
        }
    else if ("Treatment" %in% sampleSheet$condition){
        limdat.LG.CC.tw<-limdat.LG.CC[,c("ms",colnames(limdat.LG.CC)[match(sampleSheet$name[sampleSheet$condition %in% "Treatment"],colnames(limdat.LG.CC))],colnames(limdat.LG.CC)[match(sampleSheet$name[!sampleSheet$condition %in% "Treatment"],colnames(limdat.LG.CC))]),with=FALSE]
        }

    limdat.LG.CC.tw$chr<-gsub("_.+","",limdat.LG.CC.tw$ms)
    limdat.LG.CC.tw$pos<-gsub(".+_","",limdat.LG.CC.tw$ms)
    limdat.LG.CC.tw2<-limdat.LG.CC.tw[,c("chr","pos",colnames(limdat.LG.CC.tw)[2:(ncol(limdat.LG.CC.tw)-2)]),with=FALSE]
    gv<-sampleSheet$condition[match(colnames(limdat.LG.CC.tw2)[3:ncol(limdat.LG.CC.tw2)],sampleSheet$name)]###check this and modify if necessary
    ginfo<-unique(sampleSheet$condition[match(colnames(limdat.LG.CC.tw2)[3:ncol(limdat.LG.CC.tw2)],sampleSheet$name)])
    write.table(ginfo,file="groupInfo.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    cnn<-vector("numeric",length(gv))
    for(i in seq_along(ginfo)){
        cnn[which(gv %in% ginfo[i])]<-seq_along(which(gv %in% ginfo[i]))
    }

    cnv<-paste(gv,cnn,sep="_")
    colnames(limdat.LG.CC.tw2)<-c("chr","pos",cnv)

    write.table(limdat.LG.CC.tw2,file="metilene.IN.txt",sep="\t",row.names=FALSE,quote=FALSE)

}##end no CpGs passed filtering

print_sessionInfo("Analysis completed succesfully.")
