#run in R-3.3.1
Rlib<-commandArgs(trailingOnly=TRUE)[4]
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))


options(stringsAsFactors=FALSE,na.rm=TRUE)
require("limma",lib.loc=Rlib)
library("carData",lib.loc=Rlib)
require("car",lib.loc=Rlib)

###read in sample sheet

spath<-commandArgs(trailingOnly=TRUE)[2]
sampleInfo<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)


############read in methylation tracks; use (reference)-sorted bed files #############################################
#sorted bed file was used at the beginning to extract CGs from intervals of interest + 0-padding was used -> all files should have same row ordering (although possibly different from the reference)
###deal with non-0-padded tables from methylDackel and BisSNP

mpath<-commandArgs(trailingOnly=TRUE)[3]
mdir<-dir(mpath,pattern="*CpG.filt2.bed",full.names=TRUE)
mshort<-gsub(".CpG.filt2.bed","",basename(mdir))

cC<-c(rep("NULL",3),"numeric",rep("NULL",3),"character")

require(data.table,lib.loc=Rlib)

mlist<-vector("list",length(mdir))
for(i in seq_along(mdir)){
    tabi<-fread(mdir[i],colClasses=cC,sep="\t",header=TRUE)
    colnames(tabi)[colnames(tabi) %in% "Beta"]<-mshort[i]
    mlist[[i]]<-tabi
}

limdat<-Reduce(function(...) merge(..., all=T,by="ms",sort=FALSE), mlist)
limdat<-limdat[,c(1,match(sampleInfo$SampleID,colnames(limdat))),with=FALSE]


limdat.LG<-limdat 
#if(unique(grepl("MethylDackel",mdir))){
#    limdat.LG[,2:ncol(limdat.LG)]<-limdat[,2:ncol(limdat),with=FALSE]/100
#}
if("Merge" %in% colnames(sampleInfo)){
    colnames(limdat.LG)[2:ncol(limdat.LG)]<-sampleInfo$PlottingID[match(colnames(limdat.LG[,2:ncol(limdat.LG)]),sampleInfo$Merge)]}else{colnames(limdat.LG)[2:ncol(limdat.LG)]<-sampleInfo$PlottingID[match(colnames(limdat.LG[,2:ncol(limdat.LG)]),sampleInfo$SampleID)]}
limdat.LG[,2:ncol(limdat.LG)]<-limdat.LG[,2:ncol(limdat.LG)]/100
limdat.LG.CC<-limdat.LG[complete.cases(limdat.LG),] 
if(nrow(limdat.LG.CC)==0){ message("None of the single CpG sites passed the filtering.")}else{

    limdat.LG.CC[,2:ncol(limdat.LG)]<-limdat.LG.CC[,2:ncol(limdat.LG)]

    limdat.LG.CC.logit<-logit(limdat.LG.CC[,2:ncol(limdat.LG.CC),with=FALSE],percents=FALSE,adjust=0.025)
    rownames(limdat.LG.CC.logit)<-limdat.LG.CC$ms


    require("FactoMineR",lib.loc=Rlib)
    x1<-PCA(limdat.LG.CC[,-1,with=FALSE],graph=FALSE)

    pdf("limdat.LG.CC.PCA.pdf",paper="a4",bg="white")
    plot.PCA(x1,choix="var")
    dev.off()

########################## prepare density plots per group ##################################################

    require("ggplot2",lib.loc=Rlib)
    require("reshape2",lib.loc=Rlib)
    require(dplyr,lib.loc=Rlib)

#calculate and save row means
    limdat.LG.CC.L<-melt(limdat.LG.CC,id.vars="ms",value.name="Beta",variable.name="SampleID")
    limdat.LG.CC.L$SampleID<-as.character(limdat.LG.CC.L$SampleID)
    limdat.LG.CC.L$Group<-sampleInfo$Group[match(limdat.LG.CC.L$SampleID,sampleInfo$PlottingID)]
    limdat.LG.CC.Means<-data.table(summarize(group_by(limdat.LG.CC.L,ms,Group),Beta.Mean=mean(Beta)))

    head(limdat.LG.CC.Means)

##density plots
    ggplot(data=limdat.LG.CC.Means,aes(x=Beta.Mean))+geom_density(aes(group=Group,colour=Group,fill=Group),alpha=0.3)+ggtitle("Single CpG sites")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")
    ggsave(filename="Beta.MeanXgroup.all.dens.png")

##violin plots
    ggplot(data=limdat.LG.CC.Means)+geom_violin(aes(x=Group,y=Beta.Mean,fill=Group))+geom_boxplot(aes(x=Group,y=Beta.Mean),width=0.1)+ggtitle("Single CpG sites")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")
    ggsave("Beta.MeanXgroup.all.violin.png")

#limma
    design<-as.data.frame(matrix(ncol=2,nrow=(ncol(limdat.LG.CC.logit))),stringsAsFactors=FALSE)
    colnames(design)<-c("Intercept","Group")
    rownames(design)<-colnames(limdat.LG.CC.logit)
    design$Group<-as.numeric(factor(sampleInfo$Group[match(colnames(limdat.LG.CC.logit),sampleInfo$PlottingID)]))
    design$Intercept<-1
    design<-as.matrix(design)

    fit<-lmFit(limdat.LG.CC.logit,design)
    fit.eB<-eBayes(fit)
    #head(fit.eB$s2.prior)
    #head(fit.eB$s2.post)

    tT.FDR5<-topTable(fit.eB,2,p.value=0.05,number=Inf)
    if(nrow(tT.FDR5)==0){ message("No CpG sites were significantly differentially methylated.")}else{
        tT.FDR5<-tT.FDR5[,c("logFC","t","adj.P.Val","B")]
        write.table(tT.FDR5,file="limdat.LG.CC.tT.FDR5.txt",sep="\t",quote=FALSE)



        nrow(tT.FDR5)
        nrow(limdat.LG.CC.logit)
        nrow(tT.FDR5)/nrow(limdat.LG.CC.logit)

###
        limdat.LG.CC.Diff<-summarize(group_by(limdat.LG.CC.Means,ms),Diff=(Beta.Mean[1]-Beta.Mean[2]))
        tT.FDR5.Diff0.2<-tT.FDR5[rownames(tT.FDR5) %in% limdat.LG.CC.Diff$ms[abs(limdat.LG.CC.Diff$Diff)>=0.2],]

        nrow(tT.FDR5.Diff0.2)
        nrow(tT.FDR5.Diff0.2)/nrow(limdat.LG.CC.logit)

        save(limdat.LG,file="limdat.LG.RData")
        save(limdat.LG.CC.Means,limdat.LG.CC.Diff,tT.FDR5,file="singleCpG.RData")

    }###end if topTable has at least 1 entry

        limdat.LG.CC.tw<-limdat.LG.CC
        limdat.LG.CC.tw$chr<-gsub("_.+","",limdat.LG.CC.tw$ms)
        limdat.LG.CC.tw$pos<-gsub(".+_","",limdat.LG.CC.tw$ms)
        limdat.LG.CC.tw2<-limdat.LG.CC.tw[,c("chr","pos",colnames(limdat.LG.CC)[2:ncol(limdat.LG.CC)]),with=FALSE]
        gv<-sampleInfo$Group[match(colnames(limdat.LG.CC)[2:ncol(limdat.LG.CC)],sampleInfo$PlottingID)]
        gtab<-table(sampleInfo$Group[match(colnames(limdat.LG.CC)[2:ncol(limdat.LG.CC)],sampleInfo$PlottingID)])
        cnn<-vector("numeric",length(gv))
        for(i in seq_along(gtab)){
            cnn[which(gv %in% names(gtab)[i])]<-seq_along(which(gv %in% names(gtab)[i]))
        }

        cnv<-paste(gv,cnn,sep="_")
        colnames(limdat.LG.CC.tw2)<-c("chr","pos",cnv)

        write.table(limdat.LG.CC.tw2,file="metilene.IN.txt",sep="\t",row.names=FALSE,quote=FALSE)

    
}###end if any CpGs passed filtering

