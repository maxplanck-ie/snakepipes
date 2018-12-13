.libPaths(R.home("library"))

#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

options(stringsAsFactors=FALSE,na.rm=TRUE)

importfunc<-commandArgs(trailingOnly=TRUE)[9]
source(importfunc)

bedF<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("processing %s",bedF))
bedshort<-gsub(".bed","",basename(bedF))

if (length(readLines(bedF))==0) {print_sessionInfo("No DMRs found.")}else{

    bedtab<-read.table(bedF,header=FALSE,sep="\t",as.is=TRUE,quote="")
    colnames(bedtab)<-c("CHROM","START","END","qvalue","MeanDiff","NumCpGs","pMWU","p2DKS","meanA","meanB")
    bedtab$Name<-paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")
    bedtab<-unique(bedtab)

    auxF<-commandArgs(trailingOnly=TRUE)[3]
    message(sprintf("processing %s",auxF))
    auxbed<-read.table(auxF,header=FALSE,sep="\t",as.is=TRUE,quote="")
    cnpool<-c("CHROM","START","END")
    colnames(auxbed)[1:3]<-cnpool
    if(!unique(grepl("STRAND",colnames(auxbed)))){auxbed$STRAND<-"*"}
    if(!unique(grepl("Name",colnames(auxbed)))){auxbed$Name<-paste(auxbed$CHROM,auxbed$START,sep="_")}


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

    bedGR<-GRanges(seqnames=bedtab$CHROM,ranges=IRanges(start=bedtab$START,end=bedtab$END,names=paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")))
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
    if(nrow(bedtab.CC)==0) {print_sessionInfo("None of the metilene intervals passed the filtering.")}else{

        limdat.LG.inCGI<-limdat.LG.inCGI[complete.cases(limdat.LG.inCGI),] ##to be consistent with metilene input!
        CGI.limdat<-as.data.frame(apply(limdat.LG.inCGI[,2:(ncol(limdat.LG.inCGI)-3)],2,function(X){ave(X,factor(limdat.LG.inCGI$IntID),FUN=function(X)mean(X,na.rm=TRUE))}),stringsAsFactors=FALSE)

        CGI.limdat$IntID<-limdat.LG.inCGI$IntID
        CGI.limdat<-unique(CGI.limdat)
        rownames(CGI.limdat)<-CGI.limdat$IntID
        CGI.limdat<-CGI.limdat[,-grep("IntID",colnames(CGI.limdat))]

        CGI.limdat.CC<-CGI.limdat[bedtab.CC$Name,] ##this applies the bedtab$CGI.NAf filter
    ####for differential interval methylation
    ### limma + ebayes + BH p value adjustment

        require("limma")
        require("car")
        require("FactoMineR")
        require("reshape2")
        require("ggplot2")
        require("dplyr")

        CGI.limdat.CC.logit<-logit(CGI.limdat.CC,percents=FALSE,adjust=0.025)
        x1<-PCA(CGI.limdat.CC,graph=FALSE)

        if(nrow(x1$eig)>=2){
        png(paste0(bedshort,".CGI.limdat.CC.PCA.png"),bg="white") 
        plot.PCA(x1,choix="var")
        dev.off()}else{print_sessionInfo("There are not enough PC dimentions for a 2D plot.")}

    #calculate row means
        spath<-commandArgs(trailingOnly=TRUE)[5]
        sampleSheet<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)
    #calculate and save row means
        CGI.limdat.CC$IntID<-rownames(CGI.limdat.CC)
        CGI.limdat.CC.L<-melt(CGI.limdat.CC,id.vars="IntID",value.name="Beta",variable.name="SampleID")
        CGI.limdat.CC.L$Group<-sampleSheet$condition[match(CGI.limdat.CC.L$SampleID,sampleSheet$name)]
        CGI.limdat.CC.Means<-data.table(summarize(group_by(CGI.limdat.CC.L,IntID,Group),Beta.Mean=mean(Beta)))
        save(CGI.limdat.CC,file="CGI.limdat.CC.RData")

    if ("Control" %in% CGI.limdat.CC.Means$Group){
            CGI.limdat.CC.Means$Group<-factor(CGI.limdat.CC.Means$Group)
            CGI.limdat.CC.Means$Group<-relevel(CGI.limdat.CC.Means$Group,ref="Control")}
       if ("WT" %in% CGI.limdat.CC.Means$Group){
            CGI.limdat.CC.Means$Group<-factor(CGI.limdat.CC.Means$Group)
            CGI.limdat.CC.Means$Group<-relevel(CGI.limdat.CC.Means$Group,ref="WT")}

    ##density plots
        ggplot(data=CGI.limdat.CC.Means,aes(x=Beta.Mean))+geom_density(aes(group=Group,colour=Group,fill=Group),alpha=0.3)+ggtitle("Differentially methylated regions")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))+scale_colour_manual(values=c("grey28","red","darkblue","darkgreen"))+xlim(0,1)
        ggsave(paste0(bedshort,".Beta.MeanXgroup.metilene.dens.png"))

    ##violin plots
        ggplot(data=CGI.limdat.CC.Means)+geom_violin(aes(x=Group,y=Beta.Mean,fill=Group))+geom_boxplot(aes(x=Group,y=Beta.Mean),width=0.1)+ggtitle("Differentially methylated regions")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))+ylim(0,1)
        ggsave(paste0(bedshort,".Beta.MeanXgroup.metilene.violin.png"))

    #differential methylation
        design<-as.data.frame(matrix(ncol=2,nrow=(ncol(CGI.limdat.CC.logit))),stringsAsFactors=FALSE)
        colnames(design)<-c("Intercept","Group")
        rownames(design)<-colnames(CGI.limdat.CC.logit)
        if("Control" %in% sampleSheet$condition){
            gp<-factor(sampleSheet$condition[match(colnames(CGI.limdat.CC.logit),sampleSheet$name)])
            gp<-relevel(gp,ref="Control")
            design$Group<-as.numeric(gp)} else if("WT" %in% sampleSheet$condition){
            gp<-factor(sampleSheet$condition[match(colnames(CGI.limdat.CC.logit),sampleSheet$name)])
            gp<-relevel(gp,ref="WT")
            design$Group<-as.numeric(gp)}else{design$Group<-as.numeric(factor(sampleSheet$condition))}
        design$Intercept<-1
        design<-as.matrix(design)

        fit<-lmFit(CGI.limdat.CC.logit,design)
        fit.eB<-eBayes(fit)

        ##read filters from commandline args
        minAbsDiff<-as.numeric(commandArgs(trailingOnly=TRUE)[7])
        fdr<-as.numeric(commandArgs(trailingOnly=TRUE)[8])
    
        tT<-topTable(fit.eB,2,p.value=1,number=Inf)
        tT$IntID<-rownames(tT)
        plotdat<-melt(tT,measure.vars=c("P.Value","adj.P.Val"),value.name="pval",variable.name="Category",id.vars="IntID")

        ggplot(data=plotdat)+geom_histogram(aes(x=pval,group=Category,fill=Category),binwidth=0.005)+theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+scale_fill_manual(values=c("grey28","red","darkblue","darkgreen"))+geom_vline(aes(xintercept=as.numeric(fdr)))
        ggsave(paste0(bedshort,"_pvalue.distribution.png"))

### annotate top table with mean difference
        meandatW<-dcast(data=CGI.limdat.CC.Means,IntID~Group,value.var="Beta.Mean")
        if(sum(c("Control","Treatment") %in% colnames(meandatW))==2){meandatW$Diff<-with(meandatW,Treatment-Control)} else if(sum(c("WT","Mut") %in% colnames(meandatW))==2){meandatW$Diff<-with(meandatW,Mut-WT)}else{meandatW$Diff<-meandatW[,2]-meandatW[,3]}

        tT$Diff<-meandatW$Diff[match(rownames(tT),meandatW$IntID)]

        tT$Filter<-"Fail"
        tT$Filter[tT$adj.P.Val<fdr&abs(tT$Diff)>=minAbsDiff]<-"Pass"

        ggplot(data=tT)+geom_point(aes(x=Diff,y=-log10(adj.P.Val),color=Filter))+theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean difference")+scale_color_manual(values=c("grey28","red","darkblue","darkgreen"))
        ggsave(paste0(bedshort,"_volcano.plot.png"))


#### filter top table according to thresholds

        tT_filt<-tT[tT$adj.P.Val<fdr & abs(tT$Diff)>=minAbsDiff,]        

        if(nrow(tT_filt)==0){print_sessionInfo("No metilene intervals were significantly differentially methylated.")
        } else {
            tT_filt<-tT_filt[,c("logFC","t","adj.P.Val","B")]
            
            nrow(tT_filt)
            nrow(CGI.limdat.CC.logit)
            nrow(tT_filt)/nrow(CGI.limdat.CC.logit)


    #annotate metilene output with this information
            CGI.bed.intT<-as.data.frame(merge(x=bedtab,y=tT_filt,by.x="Name",by.y="row.names",sort=FALSE,all.x=TRUE))
            CGI.bed.intT<-CGI.bed.intT[,!colnames(CGI.bed.intT) %in% "Name"]
            write.table(CGI.bed.intT,file=paste0(bedshort,".limma_unfiltered.bed"),sep="\t",quote=FALSE,row.names=FALSE)
            save(CGI.bed.intT,file=paste0(bedshort,".limma_unfiltered.RData"))
            CGI.bed.intT_filt<-CGI.bed.intT[CGI.bed.intT$adj.P.Val<fdr & abs(CGI.bed.intT$MeanDiff)>=minAbsDiff,]
            CGI.bed.intT_filt<-CGI.bed.intT_filt[!is.na(CGI.bed.intT_filt$CHROM),]
            if(nrow(CGI.bed.intT_filt)>0){
            write.table(CGI.bed.intT_filt,file=paste0(bedshort,".limma_filtered.bed"),sep="\t",quote=FALSE,row.names=FALSE)}

    ####### add nearest gene information
            genMod<-commandArgs(trailingOnly=TRUE)[6]
            if (genMod!='NA' & file.exists(genMod)){
                message(sprintf("Processing gene models in %s",genMod))

                system(paste0('mkdir -p ',file.path(wdir,"temp")))

                system(paste0('bedtools sort -i ', genMod,'  > ' ,wdir ,'/temp/genes.sorted.bed'))
                system(paste0('sed -i \'/CHROM/d\' ',wdir,'/',bedshort,".limma_unfiltered.bed"))
                system(paste0('bedtools sort -i ',wdir,'/', bedshort,".limma_unfiltered.bed",' > ',wdir,'/temp/',bedshort,".limma.sorted.bed"))

                system(paste0('bedtools closest -D b -a ',wdir,'/temp/',bedshort,".limma.sorted.bed",' -b ', wdir ,'/temp/genes.sorted.bed',' > ',wdir,'/temp/',bedshort,'.limma.closest.bed'))

                DMR.filt.an<-fread(paste0(wdir,'/temp/',bedshort,'.limma.closest.bed'),header=FALSE,sep="\t")
                DMR.filt.an<-DMR.filt.an[,c(1:17,18:21,23,30),with=FALSE]
                colnames(DMR.filt.an)<-c(colnames(CGI.bed.intT),"ChrEns","StartEns","EndEns","ENST","StrandEns","Dist")

                library(biomaRt)
                emv<-c("ENSDART"="drerio","ENSMUST"="mmusculus","ENSG"="hsapiens","FBtr"="dmelanogaster")
                ems<-emv[grep(gsub("[0-9].+","",DMR.filt.an$ENST[1]),names(emv))]
                ens.xx<-useMart(biomart="ensembl",dataset=paste0(ems,"_gene_ensembl"))
                DMR.filt.an$ENST<-gsub("\\.[0-9].$","",DMR.filt.an$ENST)
                bm<-getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","external_gene_name","description"),filters="ensembl_transcript_id",values=DMR.filt.an$ENST,mart=ens.xx)

                DMR.filt.an2<-merge(x=DMR.filt.an,y=bm,by.x="ENST",by.y="ensembl_transcript_id",all.x=TRUE,allow.cartesian=TRUE)
                write.table(DMR.filt.an2,file="metilene.limma.annotated_unfiltered.txt",row.names=FALSE,quote=FALSE,sep="\t")
                DMR.filt.an2.pos<-DMR.filt.an2[DMR.filt.an2$MeanDiff>=minAbsDiff&!is.na(DMR.filt.an2$adj.P.Val),]
                if(nrow(DMR.filt.an2.pos)>0){write.table(DMR.filt.an2.pos,file="metilene.limma.annotated_filtered.UP.txt",row.names=FALSE,quote=FALSE,sep="\t")}
                DMR.filt.an2.neg<-DMR.filt.an2[DMR.filt.an2$MeanDiff<=(-minAbsDiff)&!is.na(DMR.filt.an2$adj.P.Val),] 
                if(nrow(DMR.filt.an2.neg)>0){write.table(DMR.filt.an2.neg,file="metilene.limma.annotated_filtered.DOWN.txt",row.names=FALSE,quote=FALSE,sep="\t")}
            } else {print_sessionInfo("No gene models file was provided.")}
        print_sessionInfo("Analysis completed succesfully.")
        } # end if tT_filt has at least 1 row
    } # end if any intervals passed filtering
    
} #end if bed file has at least 1 line

