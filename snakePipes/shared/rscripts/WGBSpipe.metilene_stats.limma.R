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

if (length(readLines(bedF))==0) {message("No DMRs found.")}else{

    bedtab<-read.table(bedF,header=FALSE,sep="\t",as.is=TRUE,quote="")
    colnames(bedtab)<-c("CHROM","START","END","qvalue","MeanDiff","NumCpGs","pMWU","p2DKS","meanA","meanB")
    bedtab$Name<-paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")
    bedtab<-unique(bedtab)

    auxF<-commandArgs(trailingOnly=TRUE)[3]
    message(sprintf("processing %s",auxF))
    auxbed<-read.table(auxF,header=FALSE,sep="\t",as.is=TRUE,quote="")
    cnpool<-c("CHROM","START","END")#message(sprintf("processing %s",bedF))
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
    if(nrow(bedtab.CC)==0) {message("None of the metilene intervals passed the filtering.")}else{

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
       if ("WT" %in% CGI.limdat.CC.Means$Group){
            CGI.limdat.CC.Means$Group<-factor(CGI.limdat.CC.Means$Group)
            CGI.limdat.CC.Means$Group<-relevel(CGI.limdat.CC.Means$Group,ref="WT")}

    ##density plots
        ggplot(data=CGI.limdat.CC.Means,aes(x=Beta.Mean))+geom_density(aes(group=Group,colour=Group,fill=Group),alpha=0.3)+ggtitle("Differentially methylated regions")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red"))+scale_colour_manual(values=c("grey28","red"))+xlim(0,1)
        ggsave(paste0(bedshort,".Beta.MeanXgroup.metilene.dens.png"))

    ##violin plots
        ggplot(data=CGI.limdat.CC.Means)+geom_violin(aes(x=Group,y=Beta.Mean,fill=Group))+geom_boxplot(aes(x=Group,y=Beta.Mean),width=0.1)+ggtitle("Differentially methylated regions")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")+scale_fill_manual(values=c("grey28","red"))+ylim(0,1)
        ggsave(paste0(bedshort,".Beta.MeanXgroup.metilene.violin.png"))

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
        tT.FDR5<-topTable(fit.eB,2,p.value=0.05,number=Inf)
        if(nrow(tT.FDR5)==0) {message("No metilene intervals were significantly differentially methylated.")}else{
            tT.FDR5<-tT.FDR5[,c("logFC","t","adj.P.Val","B")]
            write.table(tT.FDR5,file=paste0(bedshort,".CGI.limdat.CC.tT.FDR5.txt"),sep="\t",quote=FALSE)

            nrow(tT.FDR5)
            nrow(CGI.limdat.CC.logit)
            nrow(tT.FDR5)/nrow(CGI.limdat.CC.logit)


    #annotate metilene output with this information
            CGI.bed.intT<-as.data.frame(merge(x=bedtab,y=tT.FDR5,by.x="Name",by.y="row.names",sort=FALSE,all.x=TRUE))
            CGI.bed.intT<-CGI.bed.intT[,!colnames(CGI.bed.intT) %in% "Name"]
            write.table(CGI.bed.intT,file=paste0(bedshort,".limma.bed"),sep="\t",quote=FALSE,row.names=FALSE)
            save(CGI.bed.intT,file=paste0(bedshort,".limma.RData"))

    ####### add nearest gene information
            genMod<-commandArgs(trailingOnly=TRUE)[6]
            if (genMod!='NA' & file.exists(genMod)){
                message(sprintf("Processing gene models in %s",genMod))


                #system(paste0('sed -e \'s/^/chr/\' ', genMod,' | sort -d  -k1,1 -k2,2n  | sed -e \'s/chr//\'  > ' ,wdir ,'/genes.sorted.bed'))
                #system(paste0('sed -e \'s/^/chr/\' ',wdir,'/', bedshort,".limma.bed",' | sort  -k1,1 -k2,2n | sed -e \'s/chr//\' > ',wdir,'/',bedshort,".limma.sorted.bed"))
                #system(paste0('sed -i \'/CHROM/d\' ',wdir,'/',bedshort,".limma.sorted.bed"))

                system(paste0('bedtools sort ', genMod,'  > ' ,wdir ,'/genes.sorted.bed'))
                system(paste0('bedtools sort ',wdir,'/', bedshort,".limma.bed",' > ',wdir,'/',bedshort,".limma.sorted.bed"))
                system(paste0('sed -i \'/CHROM/d\' ',wdir,'/',bedshort,".limma.sorted.bed"))

                system(paste0('bedtools closest -D b -a ',wdir,'/',bedshort,".limma.sorted.bed",' -b ', wdir ,'/genes.sorted.bed',' > ',wdir,'/',bedshort,'.limma.closest.bed'))

                DMR.filt.an<-fread(paste0(wdir,'/',bedshort,'.limma.closest.bed'),header=FALSE,sep="\t")
                DMR.filt.an<-DMR.filt.an[,c(1:17,18:21,23,30),with=FALSE]
                colnames(DMR.filt.an)<-c(colnames(CGI.bed.intT),"ChrEns","StartEns","EndEns","ENST","StrandEns","Dist")

                library(biomaRt)
                emv<-c("ENSDART"="drerio","ENSMUST"="mmusculus","ENSG"="hsapiens","FBtr"="dmelanogaster")
                ems<-emv[grep(gsub("[0-9].+","",DMR.filt.an$ENST[1]),names(emv))]
                ens.xx<-useMart(biomart="ensembl",dataset=paste0(ems,"_gene_ensembl"))
                bm<-getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","external_gene_name","description"),filters="ensembl_transcript_id",values=DMR.filt.an$ENST,mart=ens.xx)

                DMR.filt.an2<-merge(x=DMR.filt.an,y=bm,by.x="ENST",by.y="ensembl_transcript_id",all.x=TRUE,allow.cartesian=TRUE)
                write.table(DMR.filt.an2,file="metilene.limma.annotated.txt",row.names=FALSE,quote=FALSE,sep="\t")
                DMR.filt.an2.pos<-DMR.filt.an2[DMR.filt.an2$MeanDiff>0,]
                if(nrow(DMR.filt.an2.pos)>0){write.table(DMR.filt.an2.pos,file="metilene.limma.annotated.UP.txt",row.names=FALSE,quote=FALSE,sep="\t")}
                DMR.filt.an2.neg<-DMR.filt.an2[DMR.filt.an2$MeanDiff<0,] 
                if(nrow(DMR.filt.an2.neg)>0){write.table(DMR.filt.an2.neg,file="metilene.limma.annotated.DOWN.txt",row.names=FALSE,quote=FALSE,sep="\t")}

            } else {message("No gene models file was provided.")}

        }
    }
}

