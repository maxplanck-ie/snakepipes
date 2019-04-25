.libPaths(R.home("library"))

#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

importfunc<-commandArgs(trailingOnly=TRUE)[6]
source(importfunc)


options(stringsAsFactors=FALSE,na.rm=TRUE)
require("limma")
#library("carData")
require("car")

###read in sample sheet

spath<-commandArgs(trailingOnly=TRUE)[2]
sampleSheet<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)
if(!"PlottingID" %in% colnames(sampleSheet)){sampleSheet$PlottingID<-sampleSheet$name}

require(data.table)
dpath<-commandArgs(trailingOnly=TRUE)[3]
load(file.path(dpath,"limdat.LG.RData"))
limdat.LG.CC<-limdat.LG[complete.cases(limdat.LG),] 
limdat.LG.CC.logit<-logit(limdat.LG.CC[,2:ncol(limdat.LG.CC),with=FALSE],percents=FALSE,adjust=0.025) ##result is a data.frame
rownames(limdat.LG.CC.logit)<-limdat.LG.CC$ms


require("FactoMineR")
x1<-PCA(limdat.LG.CC[,-1,with=FALSE],graph=FALSE)

png("limdat.LG.CC.PCA.png",bg="white")
plot.PCA(x1,choix="var")
dev.off()

########################## prepare density plots per group ##################################################

require("ggplot2")
require("reshape2")
require(dplyr)

#calculate and save row means
limdat.LG.CC.L<-melt(limdat.LG.CC,id.vars="ms",value.name="Beta",variable.name="SampleID")
limdat.LG.CC.L$SampleID<-as.character(limdat.LG.CC.L$SampleID)
limdat.LG.CC.L$Group<-sampleSheet$condition[match(limdat.LG.CC.L$SampleID,sampleSheet$name)]
limdat.LG.CC.Means<-data.table(summarize(group_by(limdat.LG.CC.L,ms,Group),Beta.Mean=mean(Beta)))

print(head(limdat.LG.CC.Means))

if ("Control" %in% limdat.LG.CC.Means$Group){
    limdat.LG.CC.Means$Group<-factor(limdat.LG.CC.Means$Group)
    limdat.LG.CC.Means$Group<-relevel(limdat.LG.CC.Means$Group,ref="Control")
} else if ("WT" %in% limdat.LG.CC.Means$Group){
    limdat.LG.CC.Means$Group<-factor(limdat.LG.CC.Means$Group)
    limdat.LG.CC.Means$Group<-relevel(limdat.LG.CC.Means$Group,ref="WT")
} else {
    limdat.LG.CC.Means$Group<-factor(limdat.LG.CC.Means$Group)
}

##density plots
ggplot(data=limdat.LG.CC.Means,aes(x=Beta.Mean))+geom_density(aes(group=Group,colour=Group,fill=Group),alpha=0.3)+ggtitle("Single CpG sites")+
theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))+scale_colour_manual(values=c("grey28","red","darkblue","darkgreen"))
ggsave("Beta.MeanXgroup.all.dens.png")

##violin plots
ggplot(data=limdat.LG.CC.Means)+geom_violin(aes(x=Group,y=Beta.Mean,fill=Group))+geom_boxplot(aes(x=Group,y=Beta.Mean),width=0.1)+ggtitle("Single CpG sites")+
theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14),axis.text.x = element_text(angle = 90, hjust = 1))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))
ggsave("Beta.MeanXgroup.all.violin.png")

print(levels(limdat.LG.CC.Means$Group))
if(length(levels(limdat.LG.CC.Means$Group))==2){

    sink("GroupMean.ttest.txt")
    limdat.LG.CC.MeansXSample<-data.table(summarize(group_by(limdat.LG.CC.L,SampleID),Beta.Mean=mean(Beta)))
    limdat.LG.CC.MeansXSample$Beta.Mean.logit<-logit(limdat.LG.CC.MeansXSample$Beta.Mean,percents=FALSE,adjust=0.025)
    limdat.LG.CC.MeansXSample$Group<-sampleSheet$condition[match(limdat.LG.CC.MeansXSample$SampleID,sampleSheet$name)]
    print(limdat.LG.CC.MeansXSample)
    print(t.test(Beta.Mean.logit~Group,data=limdat.LG.CC.MeansXSample,var.equal=TRUE))
    sink()



    #limma
    design<-as.data.frame(matrix(ncol=2,nrow=(ncol(limdat.LG.CC.logit))),stringsAsFactors=FALSE)
    colnames(design)<-c("Intercept","Group")
    rownames(design)<-colnames(limdat.LG.CC.logit)
    if("Control" %in% sampleSheet$condition){
        gp<-factor(sampleSheet$condition[match(colnames(limdat.LG.CC.logit),sampleSheet$name)])
        gp<-relevel(gp,ref="Control")
        design$Group<-as.numeric(gp)}
    else if("WT" %in% sampleSheet$condition){
        gp<-factor(sampleSheet$condition[match(colnames(limdat.LG.CC.logit),sampleSheet$name)])
        gp<-relevel(gp,ref="WT")
        design$Group<-as.numeric(gp)}
    else{design$Group<-as.numeric(factor(sampleSheet$condition))}
    design$Intercept<-1
    design<-as.matrix(design)

    fit<-lmFit(limdat.LG.CC.logit,design)
    fit.eB<-eBayes(fit)

##read filters from commandline args
    minAbsDiff<-commandArgs(trailingOnly=TRUE)[4]
    fdr<-commandArgs(trailingOnly=TRUE)[5]
    
    tT<-topTable(fit.eB,2,p.value=1,number=Inf)
    tT$IntID<-rownames(tT)
    plotdat<-melt(tT,measure.vars=c("P.Value","adj.P.Val"),value.name="pval",variable.name="Category",id.vars="IntID")

    ggplot(data=plotdat)+geom_histogram(aes(x=pval,group=Category,fill=Category),binwidth=0.005)+theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))+geom_vline(aes(xintercept=as.numeric(fdr)))
    ggsave("SingleCpG_pvalue.distribution.png")

### annotate top table with mean difference
    meandatW<-dcast(data=limdat.LG.CC.Means,ms~Group,value.var="Beta.Mean")
    if(sum(c("Control","Treatment") %in% colnames(meandatW))==2){meandatW$Diff<-with(meandatW,Treatment-Control)}
    if(sum(c("WT","Mut") %in% colnames(meandatW))==2){meandatW$Diff<-with(meandatW,Mut-WT)}else{meandatW$Diff<-meandatW[,2]-meandatW[,3]}
    tT$Diff<-meandatW$Diff[match(tT$IntID,meandatW$ms)]

    tT$Filter<-"Fail"
    tT$Filter[tT$adj.P.Val<fdr&abs(tT$Diff)>=minAbsDiff]<-"Pass"

    ggplot(data=tT)+geom_point(aes(x=Diff,y=-log10(adj.P.Val),color=Filter))+theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean difference")+scale_color_manual(values=c("grey28","red","darkblue","darkgreen"))
    ggsave("SingleCpG_volcano.plot.png")


#### filter top table according to thresholds

    tT_filt<-tT[tT$adj.P.Val<fdr & abs(tT$Diff)>=minAbsDiff,] 

    if(nrow(tT_filt)==0){ print_sessionInfo("No CpG sites were significantly differentially methylated according to the thresholds.")}else{
            
        nrow(tT_filt)
        nrow(limdat.LG.CC.logit)
        nrow(tT_filt)/nrow(limdat.LG.CC.logit)

        save(limdat.LG.CC.Means,tT_filt,file="singleCpG.RData")
        print_sessionInfo("Analysis completed succesfully.")

    }###end if topTable has at least 1 entry


} else {print_sessionInfo('More than 2 sample groups were provided. No statistical inference will be computed.')}### end if exactly two sample groups were specified

  


