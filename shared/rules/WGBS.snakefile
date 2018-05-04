import os
import re
from operator import is_not

###get automatic cut threshold for hard-trimming of 5' ends

if trimReads=='auto':
    if fqcin:
        rule get_cut_thd:
            input:
                R1zip = fqcin+"/{sample}"+reads[0]+"_fastqc.zip",
                R2zip = fqcin+"/{sample}"+reads[1]+"_fastqc.zip"
            output:
                R12ct= "FastQC_In/{sample}"+".R12.ct.txt"
            log:"FastQC_In/logs/{sample}.cutThd.log"
            run:
                for f,g,z,l in zip ({input.R1zip},{input.R2zip},output,log):
                    with open(z, "w") as oo:
                        cutThdRes=calc_cutThd([f,g],fqcin,l,wdir)
                        oo.write('\n'.join('%s\t%s\n' % x for x in cutThdRes))
                os.chdir(wdir)


    rule trimReads:
        input:
            R1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/{sample}"+reads[1]+".fastq.gz",
            R12ct= "FastQC_In/{sample}"+".R12.ct.txt"
        output:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        #log:"FASTQ_Cutadapt/logs/{sample}.log"
        run:
            shell(cut_reads_auto(input.R12ct,nthreads))

elif trimReads=='user':
    rule trimReads:
        input:
            R1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        #log:"FASTQ_Cutadapt/logs/{sample}.log"
        run:
            shell(cut_reads_user(trimThreshold,nthreads,trimOtherArgs,nextera))



if not trimReads is None:
    rule postTrimFastQC:
        input:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        output:
            R1fqc="FastQC_Cutadapt/{sample}"+reads[0]+"_fastqc.html",
            R2fqc="FastQC_Cutadapt/{sample}"+reads[1]+"_fastqc.html"
        #log:"FastQC_Cutadapt/logs/{sample}.log"
        params:
            fqcout=os.path.join(wdir,'FastQC_Cutadapt')
        run:
            shell(os.path.join(fastqc_path,'fastqc ')+' --outdir ' + "{params.fqcout}" + ' -t  ' + str(nthreads) + ' ' + "{input.R1cut}" + ' ' + "{input.R2cut}")

if convRef:
    rule conv_ref:
        input:
            refG=refG
        output:
            crefG=os.path.join('conv_ref',re.sub('.fa','.fa.bwameth.c2t',os.path.basename(refG)))
        #log:"conv_ref/logs/convref.log"
        shell:'ln -s '+ "{input.refG}" + ' '+ os.path.join("aux_files",os.path.basename("{input.refG}"))+'; python ' + os.path.join(bwameth_path,'bwameth.py') + ' index ' + os.path.join("aux_files",os.path.basename("{input.refG}"))

if not trimReads is None:
    rule map_reads:
        input:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz",
            crefG=crefG
        output:
            temp(sbam="bams/{sample}"+".sorted.bam")
        #log:"bams/logs/{sample}"+".mapping.log"
        run:
            shell(bMeth_map_reads(input.R1cut,input.R2cut,input.crefG,nthreads))

if trimReads is None:
    rule map_reads:
        input:
            R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2="FASTQ/{sample}"+reads[1]+".fastq.gz",
            crefG=crefG
        output:
            sbam="bams/{sample}"+".sorted.bam"
        #log:"bams/logs/{sample}"+".mapping.log"
        run:
            shell(bMeth_map_reads(input.R1,input.R2,input.crefG,nthreads))


rule index_bam:
    input:
        sbam="bams/{sample}"+".sorted.bam"
    output:
        temp(sbami="bams/{sample}"+".sorted.bam.bai")
    params:
    #log:"bams/logs/{sample}"+".indexing.log"
    run:
        shell(os.path.join(samtools_path,'samtools') + ' index ' + "{input.sbam}")

rule rm_dupes:
    input:
        sbami="bams/{sample}"+".sorted.bam.bai",
        sbam="bams/{sample}"+".sorted.bam"
    output:
        rmDupbam="bams/{sample}"+".PCRrm.bam"
    #log:"bams/logs/{sample}"+".PCRdupRm.log"
    run:
        shell(sambamba_module_path + 'sambamba_v0.6.6 markdup --remove-duplicates -t ' + str(nthreads) + ' --tmpdir ' + tempfile.mkdtemp(suffix='',prefix='',dir='/data/extended') +' ' + "{input.sbam}" + ' ' + "{output.rmDupbam}") 

rule index_PCRrm_bam:
    input:
        sbam="bams/{sample}"+".PCRrm.bam"
    output:
        sbami="bams/{sample}"+".PCRrm.bam.bai"
    params:
    #log:"bams/logs/{sample}"+".indexing.log"
    run:
        shell(os.path.join(samtools_path,'samtools') + ' index ' + "{input.sbam}")


rule get_ran_CG:
    input:
        refG=refG
    output:
        pozF="aux_files/"+re.sub('.fa*','.poz.gz',os.path.basename(refG)),
        ranCG=os.path.join("aux_files",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
    #log:"aux_files/genome.get_ranCG.log"
    run:
        shell(mCT_get_ranCG())


rule calc_Mbias:
    input:
        refG=refG,
        rmDupBam="bams/{sample}"+".PCRrm.bam",
        sbami="bams/{sample}"+".PCRrm.bam.bai"
    output:
        mbiasTXT="QC_metrics/{sample}"+".Mbias.txt"
    #log:"QC_metrics/logs/{sample}"+".mbias.log"
    run:
        shell(os.path.join(POM_path,'MethylDackel') + ' mbias ' + "{input.refG}" + ' ' + "{input.rmDupBam}" +' ' + "{output.mbiasTXT}" +' -@ '+str(nthreads) +' 2>' + "{output.mbiasTXT}")


if intList:
    rule depth_of_cov:
        input:
            refG=refG,
            rmDupBam="bams/{sample}"+".PCRrm.bam",
            sbami="bams/{sample}"+".PCRrm.bam.bai"
        output:
            outFileList=calc_doc(intList,True)
        params:
            auxdir=os.path.join(wdir,"aux_files")
        #log:"QC_metrics/logs/{sample}"+".doc.log"
        run:
            shell(BS_doc_XT(output.outFileList,intList,input.refG))

else:
    rule depth_of_cov:
        input:
            refG=refG,
            rmDupBam="bams/{sample}"+".PCRrm.bam",
            sbami="bams/{sample}"+".PCRrm.bam.bai"
        output:
            outFileList=calc_doc(intList,True)
        params:
            auxdir=os.path.join(wdir,"aux_files")
        #log:"QC_metrics/logs/{sample}"+".doc.log"
        run:
            shell(BS_doc(output.outFileList,input.refG))

if not trimReads is None:
    rule conv_rate:
        input:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        output:
            R12cr="QC_metrics/{sample}" + '.conv.rate.txt'
        #log:"QC_metrics/logs/{sample}"+".convrate.log"
        run:
            shell(os.path.join(workflow_tools,'conversionRate_KS.sh ')+ os.path.join('FASTQ_Cutadapt',re.sub('_R1.fastq.gz','',os.path.basename(input.R1cut))) + ' ' + "{output.R12cr}")

else:
    rule conv_rate:
        input:
            R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2="FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R12cr="QC_metrics/{sample}" + '.conv.rate.txt'
        #log:"QC_metrics/logs/{sample}"+".convrate.log"
        run:
            shell(os.path.join(workflow_tools,'conversionRate_KS.sh ')+ os.path.join('FASTQ',re.sub('_R1.fastq.gz','',os.path.basename(input.R1))) + ' ' + "{output.R12cr}")

rule get_flagstat:
    input:
        rmDupbam="bams/{sample}"+".PCRrm.bam"
    output:
        fstat="QC_metrics/{sample}" + '.flagstat'
    #log:"QC_metrics/logs/{sample}"+".flagstat.log"
    run:
        shell(os.path.join(samtools_path,'samtools') + ' flagstat ' + "{input.rmDupbam}" +' > ' + "{output.fstat}") 

rule produce_report:
    input:
        doc_res=calc_doc(intList,False),
        R12cr=expand("QC_metrics/{sample}.conv.rate.txt",sample=samples),
        mbiasTXT=expand("QC_metrics/{sample}.Mbias.txt",sample=samples),
        fstat=expand("QC_metrics/{sample}.flagstat",sample=samples)
    output:
        QCrep='QC_metrics/QC_report.pdf'
    params:
        auxdir=os.path.join(wdir,"aux_files")
    #log:"QC_metrics/logs/QCrep.log"
    run:
        shell(BS_QC_rep())  

rule methyl_extract:
    input:
        rmDupbam="bams/{sample}"+".PCRrm.bam",
        sbami="bams/{sample}"+".PCRrm.bam.bai",
        refG=refG,
        mbiasTXT="QC_metrics/{sample}"+".Mbias.txt"
    output:
        methTab="methXT/{sample}_CpG.bedGraph"
    params:
        mbias_action=mbias_ignore,
        auxdir=os.path.join(wdir,"aux_files")
    #log:"methXT/logs/{sample}.methXT.log"
    run:
        shell(methXT_POM(input.rmDupbam,params.mbias_action,nthreads,wdir,output.methTab))

rule CpG_filt:
    input:
        methTab="methXT/{sample}_CpG.bedGraph"
    output:
        tabFilt="methXT/{sample}.CpG.filt2.bed"
    params:
        Rlib=R_libs_path,
        methDir=os.path.join(wdir,"methXT"),
        blacklistF=blackList
    #log:"methXT/logs/{sample}.CpG.filt.log"
    run:
        shell(filt_POM(input.methTab,params.blacklistF))


if sampleInfo:
    rule CpG_stats:
        input: expand("methXT/{sample}.CpG.filt2.bed",sample=samples)
        output:
            RDatAll='singleCpG_stats_limma/singleCpG.RData',
            Limdat='singleCpG_stats_limma/limdat.LG.RData',
            MetIN='singleCpG_stats_limma/metilene.IN.txt'
        params:
            statdir=os.path.join(wdir,'singleCpG_stats_limma'),
            sampleInfo=sampleInfo,
            Rlib=R_libs_path
        #log:"singleCpG_stats_limma/logs/singleCpGstats.log"
        run:
            shell(os.path.join(R_path,'Rscript') +' --no-save --no-restore ' + os.path.join(workflow_tools,'WGBSpipe.singleCpGstats.limma.R ') + "{params.statdir}" + ' ' + "{params.sampleInfo}" + ' ' + os.path.join(wdir,"methXT") + ' ' + "{params.Rlib}" +' ;sleep 300')

    rule run_metilene:
        input:
            MetIN='singleCpG_stats_limma/metilene.IN.txt',
            sampleInfo=sampleInfo
        output:
            MetBed='metilene_out/singleCpG.metilene.bed'
        params:
            DMRout=os.path.join(wdir,'metilene_out')
        #log:"metilene_out/logs/metilene.log"
        run:
            shell(DMR_metilene(nthreads,input.sampleInfo))

    rule get_CG_metilene:
        input:
            refG=refG,
            MetBed='metilene_out/singleCpG.metilene.bed',
            pozF="aux_files/"+re.sub('.fa*','.poz.gz',os.path.basename(refG))
        output:
            MetCG=os.path.join("aux_files",re.sub('.fa','.metilene.CpGlist.bed',os.path.basename(refG)))
        params:
            auxdir=os.path.join(wdir,"aux_files")
        #log:"aux_files/getCG_metilene.log"
        run:
            shell(mCT_get_CpGxInt(input.pozF,["{output.MetCG}"],["{input.MetBed}"],wdir))

    rule cleanup_metilene:
        input:
            Limdat='singleCpG_stats_limma/limdat.LG.RData',
            MetBed='metilene_out/singleCpG.metilene.bed',
            MetCG=os.path.join("aux_files",re.sub('.fa','.metilene.CpGlist.bed',os.path.basename(refG))),
            sampleInfo=sampleInfo,
            refG=refG
        output:
            LimBed='metilene_out/singleCpG.metilene.limma.bed',
            LimAnnot='metilene_out/metilene.limma.annotated.txt'
        params:
            Rlib=R_libs_path,
            DMRout=os.path.join(wdir,'metilene_out')
        #log:'metilene_out/logs/cleanup.log'
        run:
            shell(clean_up_metilene())


if intList:
    rule get_CG_per_int:
        input:
            intList=intList,
            refG=refG,
            pozF="aux_files/"+re.sub('.fa*','.poz.gz',os.path.basename(refG))
        output:
            outList=run_int_aggStats(intList,False) #check for full path
        #log:"aux_files/getCG_intList.log"
        run:
            shell(mCT_get_CpGxInt(input.pozF,[x for x in output.outList],[y for y in input.intList],wdir))

    if sampleInfo:
        rule intAgg_stats:
            input:
                Limdat='singleCpG_stats_limma/limdat.LG.RData',
                intList=intList,
                refG=refG,
                sampleInfo=sampleInfo
            output:
                outFiles=run_int_aggStats(intList,sampleInfo)
            params:
                auxdir=os.path.join(wdir,"aux_files"),
                aggStatdir=os.path.join(wdir,'aggregate_stats_limma'),
                Rlib=R_libs_path
            #log:'aggregate_stats_limma/logs/aggStats.limma.log'
            run:
                shell(int_stats_limma([y for y in input.intList]))

 
