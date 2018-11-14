#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

options(stringsAsFactors=FALSE,na.rm=TRUE)

bedF<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("processing %s",bedF))
bedshort<-gsub(".bed","",basename(bedF))

bedtab<-read.table(bedF,header=FALSE,sep="\t",as.is=TRUE,quote="")
cnpool<-c("CHROM","START","END")#message(sprintf("processing %s",bedF))
bedshort<-gsub(".bed","",basename(bedF))
bedtab<-read.table(bedF,header=FALSE,sep="\t",as.is=TRUE,quote="")
colnames(bedtab)[1:3]<-cnpool#[1:ncol(bedtab)]
if(!unique(grepl("STRAND",colnames(bedtab)))){bedtab$STRAND<-"*"}
if(!unique(grepl("Name",colnames(bedtab)))){bedtab$Name<-paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")}
bedtab<-unique(bedtab)

auxF<-commandArgs(trailingOnly=TRUE)[3]
message(sprintf("processing %s",auxF))
auxbed<-read.table(auxF,header=FALSE,sep="\t",as.is=TRUE,quote="")
colnames(auxbed)[1:3]<-cnpool#[1:ncol(bedtab)]
if(!unique(grepl("STRAND",colnames(auxbed)))){auxbed$STRAND<-"*"}
if(!unique(grepl("Name",colnames(auxbed)))){auxbed$Name<-paste(auxbed$CHROM,auxbed$START,sep="_")}#auxbed$END

#load single CpG data, count NAs per row/CpG
sCpGF<-commandArgs(trailingOnly=TRUE)[4]
message(sprintf("loading %s",sCpGF))

require(data.table)
load(sCpGF)

noNA<-apply(limdat.LG,1,function(X)sum(is.na(X)))
NAf<-ifelse(noNA>1,1,0)
limdat.LG$NAf<-NAf
colnames(limdat.LG)[colnames(limdat.LG) %in% "ms"]<-"Name"

require(GenomicRanges)

bedGR<-GRanges(seqnames=bedtab$CHROM,ranges=IRanges(start=bedtab$START,end=bedtab$END,names=bedtab$Name),strand=bedtab$STRAND)
auxdat<-as.data.table(merge(x=auxbed,y=limdat.LG,by="Name",all.x=TRUE,sort=FALSE))
auxdat<-auxdat[,!colnames(auxdat) %in% c("CHROM","START","END","V4","V5","STRAND"),with=FALSE]

auxGR<-GRanges(seqnames=auxbed$CHROM,ranges=IRanges(start=auxbed$START,end=auxbed$END,names=auxbed$Name),strand=auxbed$STRAND)
auxdat<-auxdat[match(names(auxGR),auxdat$Name),] ##there is some reordering during GRanges construction
############################# UPDATED CODE
auxmtch<-as.data.frame(findOverlaps(query=auxGR, subject=bedGR))

limdat.LG.inCGI<-auxdat[auxmtch$queryHits,] #for compatibility with previous method; not really necessary as auxdat comes from a pre-intersected bed file
limdat.LG.inCGI$IntID<-names(ranges(bedGR))[auxmtch$subjectHits] ##needed for grouping per interval
##
noNA2<-apply(limdat.LG.inCGI[,!colnames(limdat.LG.inCGI) %in% c("Name","NAf","IntID"),with=FALSE],1,function(X)sum(is.na(X)))
NAf2<-ifelse(noNA2>1,1,0)
limdat.LG.inCGI$NAf<-NAf2

########################### END UPDATED CODE
NA.inCGI<-with(limdat.LG.inCGI,ave(NAf,factor(IntID),FUN=sum))
limdat.LG.inCGI$NA.inCGI<-NA.inCGI

to.m<-limdat.LG.inCGI[,c("IntID", "NA.inCGI"),with=FALSE]

CGI.map<-unique(to.m)
bedtab$N.CG.NA<-CGI.map$NA.inCGI[match(bedtab$Name,CGI.map$IntID)]

N.CG.tot<-with(limdat.LG.inCGI,ave(NAf,IntID,FUN=length))
bedtab$N.CG.tot<-N.CG.tot[match(bedtab$Name,limdat.LG.inCGI$IntID)]

bedtab$CGI.NAf<-ifelse(bedtab$N.CG.NA>(0.8*bedtab$N.CG.tot),NA,1)
bedtab.CC<-bedtab[complete.cases(bedtab),]
if(nrow(bedtab.CC)==0) {message("None of the genomic intervals passed the filtering.")}else{

    CGI.limdat<-as.data.frame(apply(limdat.LG.inCGI[,2:(ncol(limdat.LG.inCGI)-3)],2,function(X){ave(X,factor(limdat.LG.inCGI$IntID),FUN=function(X)mean(X,na.rm=TRUE))}),stringsAsFactors=FALSE)

    CGI.limdat$IntID<-limdat.LG.inCGI$IntID
    CGI.limdat<-unique(CGI.limdat)
    rownames(CGI.limdat)<-CGI.limdat$IntID
    CGI.limdat<-CGI.limdat[,-grep("IntID",colnames(CGI.limdat))]

    CGI.limdat.CC<-CGI.limdat[bedtab.CC$Name,]


####for differential interval methylation
### limma + ebayes + BH p value adjustment

    require("limma")
    #library("carData")
    require("car")
    require("FactoMineR")
    require("reshape2")
    require("ggplot2")
    require("dplyr")

    CGI.limdat.CC.logit<-logit(CGI.limdat.CC,percents=FALSE,adjust=0.025)
    x1<-PCA(CGI.limdat.CC,graph=FALSE)

    if(nrow(x1$eig)>=2){
    pdf(paste0(bedshort,".CGI.limdat.CC.PCA.pdf"),paper="a4",bg="white") 
    plot.PCA(x1,choix="var")
    dev.off()}else{message("There are not enough PC dimentions for a 2D plot.")}

#calculate row means
    spath<-commandArgs(trailingOnly=TRUE)[5]
    sampleInfo<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)
#calculate and save row means
    CGI.limdat.CC$IntID<-rownames(CGI.limdat.CC)
    CGI.limdat.CC.L<-melt(CGI.limdat.CC,id.vars="IntID",value.name="Beta",variable.name="SampleID")
    CGI.limdat.CC.L$Group<-sampleInfo$Group[match(CGI.limdat.CC.L$SampleID,sampleInfo$SampleID)]
    CGI.limdat.CC.Means<-data.table(summarize(group_by(CGI.limdat.CC.L,IntID,Group),Beta.Mean=mean(Beta)))


   if ("Control" %in% CGI.limdat.CC.Means$Group){
        CGI.limdat.CC.Means$Group<-factor(CGI.limdat.CC.Means$Group)
        CGI.limdat.CC.Means$Group<-relevel(CGI.limdat.CC.Means$Group,ref="Control")}
   else if ("WT" %in% CGI.limdat.CC.Means$Group){
        CGI.limdat.CC.Means$Group<-factor(CGI.limdat.CC.Means$Group)
        CGI.limdat.CC.Means$Group<-relevel(CGI.limdat.CC.Means$Group,ref="WT")}
   else {CGI.limdat.CC.Means$Group<-factor(CGI.limdat.CC.Means$Group)}


##density plots
    ggplot(data=CGI.limdat.CC.Means,aes(x=Beta.Mean))+geom_density(aes(group=Group,colour=Group,fill=Group),alpha=0.3)+ggtitle("Genomic intervals")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))+scale_colour_manual(values=c("grey28","red","darkblue","darkgreen"))
    ggsave(paste0(bedshort,".Beta.MeanXgroup.int.dens.png"))

##violin plots
    ggplot(data=CGI.limdat.CC.Means)+geom_violin(aes(x=Group,y=Beta.Mean,fill=Group))+geom_boxplot(aes(x=Group,y=Beta.Mean),width=0.1)+ggtitle("Genomic intervals")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14),axis.text.x = element_text(angle = 90, hjust = 1))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))
    ggsave(paste0(bedshort,".Beta.MeanXgroup.int.violin.png"))

    if(length(levels(CGI.limdat.CC.Means$Group))==2){
#differential methylation
        design<-as.data.frame(matrix(ncol=2,nrow=(ncol(CGI.limdat.CC.logit))),stringsAsFactors=FALSE)
        colnames(design)<-c("Intercept","Group")
        rownames(design)<-colnames(CGI.limdat.CC.logit)
        if("Control" %in% sampleInfo$Group){
            gp<-factor(sampleInfo$Group[match(colnames(CGI.limdat.CC.logit),sampleInfo$SampleID)])
            gp<-relevel(gp,ref="Control")
            design$Group<-as.numeric(gp)}
        if("WT" %in% sampleInfo$Group){
            gp<-factor(sampleInfo$Group[match(colnames(CGI.limdat.CC.logit),sampleInfo$SampleID)])
            gp<-relevel(gp,ref="WT")
            design$Group<-as.numeric(gp)}
        else{design$Group<-as.numeric(factor(sampleInfo$Group))}
        design$Intercept<-1
        design<-as.matrix(design)

        fit<-lmFit(CGI.limdat.CC.logit,design)
        fit.eB<-eBayes(fit)


        ##read filters from commandline args
        minAbsDiff<-commandArgs(trailingOnly=TRUE)[6]
        fdr<-commandArgs(trailingOnly=TRUE)[7]
    
        tT<-topTable(fit.eB,2,p.value=1,number=Inf)
        tT$IntID<-rownames(tT)
        plotdat<-melt(tT,measure.vars=c("P.Value","adj.P.Val"),value.name="pval",variable.name="Category",id.vars="IntID")

        ggplot(data=plotdat)+geom_histogram(aes(x=pval,group=Category,fill=Category),binwidth=0.005)+theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))+geom_vline(aes(xintercept=0.02))
        ggsave(paste0(bedshort,"_pvalue.distribution.png"))

### annotate top table with mean difference
        meandatW<-dcast(data=CGI.limdat.CC.Means,ms~Group,value.var="Beta.Mean")
        if(sum(c("Control","Treatment") %in% colnames(meandatW))==2){meandatW$Diff<-with(meandatW,Treatment-Control)}
        if(sum(c("WT","Mut") %in% colnames(meandatW))==2){meandatW$Diff<-with(meandatW,Mut-WT)}
        else{meandatW$Diff<-meandatW[2]-meandatW[3]}

        tT$Diff<-meandatW$Diff[match(rownames(tT),meandatW$ms)]

        tT$Filter<-"Fail"
        tT$Filter[tT$adj.P.Val<fdr&abs(tT$Diff)>=minAbsDiff]<-"Pass"

        ggplot(data=tT)+geom_point(aes(x=Diff,y=-log10(adj.P.Val),color=Filter))+theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean difference")+scale_color_manual(values=c("grey28","red","darkblue","darkgreen"))
        ggsave(paste0(bedshort,"_volcano.plot.png"))


#### filter top table according to thresholds

        tT_filt<-tT[tT$adj.P.Val<fdr & abs(tT$Diff)>=minAbsDiff,]  

        if(nrow(tT_filt)==0) {message("No genomic intervals were significantly differentially methylated.")
            save(bedtab,limdat.LG.inCGI,CGI.limdat.CC,CGI.limdat.CC.Means,file=paste0(bedshort,".aggCpG.RData"))

        }else{
            write.table(tT_filt,file=paste0(bedshort,".CGI.limdat.CC.tT_filt.txt"),sep="\t",quote=FALSE)

            nrow(tT_filt)
            nrow(CGI.limdat.CC.logit)
            nrow(tT_filt)/nrow(CGI.limdat.CC.logit)

            save(bedtab,limdat.LG.inCGI,CGI.limdat.CC,CGI.limdat.CC.Means,tT_filt,file=paste0(bedshort,".aggCpG.RData"))

        }###end if topTable has at least 1 entry
    } else {message('More than 2 sample groups were provided. No statistical inference will be computed.')}### end if exactly two sample groups were specified

}

