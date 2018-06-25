import os
import re
from operator import is_not
import tempfile
import pandas

###conda environments:
CondaEnvironment='envs/WGBSconda.yml'
mCtCondaEnvironment='envs/methylCtools.yml'
RmdCondaEnvironment='envs/Rmdconda.yml'
CONDA_SHARED_ENV = "envs/shared_environment.yaml"

  
###get automatic cut threshold for hard-trimming of 5' ends
if trimReads=='auto':
    if fqcin:
        rule get_cut_thd:
            input:
                R1zip = fqcin+"/{sample}"+reads[0]+"_fastqc.zip",
                R2zip = fqcin+"/{sample}"+reads[1]+"_fastqc.zip"
            output:
                R12ct= "FastQC_In/{sample}.R12.ct.txt"
            log:"FastQC_In/logs/{sample}.get_cut_thd.log"
            threads: 1
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
            R12ct= "FastQC_In/{sample}.R12.ct.txt"
        output:
            R1cut=temp("FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz"),
            R2cut=temp("FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz")
        log:
            err="FASTQ_Cutadapt/logs/{sample}.trimReads.err",
            out="FASTQ_Cutadapt/logs/{sample}.trimReads.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: "ct=($(cat {input.R12ct})); cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5 -j {threads} -u ${{ct[0]}}  -U ${{ct[1]}} -o {output.R1cut} -p {output.R2cut} {input.R1} {input.R2} 1>{log.out} 2>{log.err}"
            

elif trimReads=='user':
    rule trimReads:
        input:
            R1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R1cut=temp("FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz"),
            R2cut=temp("FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz")
        log:
            err="FASTQ_Cutadapt/logs/{sample}.trimReads.err",
            out="FASTQ_Cutadapt/logs/{sample}.trimReads.out"
        params:
            adapterSeq=adapterSeq,
            trimThreshold=trimThreshold,
            trimOtherArgs=trimOtherArgs
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: "cutadapt -a {params.adapterSeq} -A {params.adapterSeq} -q {params.trimThreshold} -m 30 -j {threads} {params.trimOtherArgs} -o {output.R1cut} -p {output.R2cut} {input.R1} {input.R2} 1>{log.out} 2>{log.err}" 



if not trimReads is None:
    rule postTrimFastQC:
        input:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        output:
            R1fqc="FastQC_Cutadapt/{sample}"+reads[0]+"_fastqc.html",
            R2fqc="FastQC_Cutadapt/{sample}"+reads[1]+"_fastqc.html"
        log:
            err="FastQC_Cutadapt/logs/{sample}.postTrimFastQC.err",
            out="FastQC_Cutadapt/logs/{sample}.postTrimFastQC.out"
        params:
            fqcout=os.path.join(wdir,'FastQC_Cutadapt')
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: "fastqc --outdir {params.fqcout} -t  {threads} {input.R1cut} {input.R2cut} 1>{log.out} 2>{log.err}"

if convRef:
    rule conv_ref:
        input:
            refG=refG
        output:
            crefG=temp(os.path.join('conv_ref',re.sub('.fa','.fa.bwameth.c2t',os.path.basename(refG))))
        log:
            err="conv_ref/logs/{sample}.conv_ref.err",
            out="conv_ref/logs/{sample}.conv_ref.out"
        threads: 1
        conda: CondaEnvironment
        shell:"ln -s {input.refG} " + os.path.join("aux_files",os.path.basename("{input.refG}"))+';bwameth index ' + os.path.join("aux_files",os.path.basename("{input.refG}")) + " 1>{log.out} 2>{log.err}"

if not trimReads is None:
    rule map_reads:
        input:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz",
            crefG=crefG
        output:
            sbam=temp("bams/{sample}.sorted.bam")
        log:
            err="bams/logs/{sample}.map_reads.err",
            out="bams/logs/{sample}.map_reads.out"
        params:
            tempdir=tempfile.mkdtemp(suffix='',prefix="{sample}",dir=tempdir),
            sortThreads=min(nthreads,4),
            RG=lambda wildcards: RG_dict[wildcards.sample]
        threads: nthreads
        conda: CondaEnvironment
        shell: "bwameth.py --threads  {threads}  --read-group {params.RG} --reference {input.crefG} {input.R1cut} {input.R2cut} | samtools sort -T {params.tempdir} -m 3G -@ {params.sortThreads} -o {output.sbam} 1>{log.out} 2>{log.err}"

if trimReads is None:
    rule map_reads:
        input:
            R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2="FASTQ/{sample}"+reads[1]+".fastq.gz",
            crefG=crefG
        output:
            sbam=temp("bams/{sample}.sorted.bam")
        log:
            err="bams/logs/{sample}.map_reads.err",
            out="bams/logs/{sample}.map_reads.out"
        params:
            tempdir=tempfile.mkdtemp(suffix='',prefix="{sample}",dir=tempdir),
            sortThreads=min(nthreads,4),
            RG=lambda wildcards: RG_dict[wildcards.sample]
        threads: nthreads
        conda: CondaEnvironment
        shell: "bwameth.py --threads  {threads}  --read-group {params.RG} --reference {input.crefG} {input.R1} {input.R2} | samtools sort -T {params.tempdir} -m 3G -@ {params.sortThreads} -o {output.sbam} 1>{log.out} 2>{log.err}"


rule index_bam:
    input:
        sbam="bams/{sample}.sorted.bam"
    output:
        sbami=temp("bams/{sample}.sorted.bam.bai")
    log:
        err="bams/logs/{sample}.index_bam.err",
        out="bams/logs/{sample}.index_bam.out"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input.sbam} 1>{log.out} 2>{log.err}"

rule rm_dupes:
    input:
        sbami="bams/{sample}.sorted.bam.bai",
        sbam="bams/{sample}.sorted.bam"
    output:
        rmDupbam="bams/{sample}.PCRrm.bam"
    log:
        err="bams/logs/{sample}.rm_dupes.err",
        out="bams/logs/{sample}.rm_dupes.out"
    params:
        tempdir=tempfile.mkdtemp(suffix='',prefix='',dir=tempdir)
    threads: nthreads
    conda: CONDA_SHARED_ENV
    shell: "sambamba markdup --remove-duplicates -t {threads} --tmpdir {params.tempdir} {input.sbam} {output.rmDupbam} 1>{log.out} 2>{log.err}" 

rule index_PCRrm_bam:
    input:
        sbam="bams/{sample}.PCRrm.bam"
    output:
        sbami="bams/{sample}.PCRrm.bam.bai"
    params:
    log:
        err="bams/logs/{sample}.index_PCRrm_bam.err",
        out="bams/logs/{sample}.index_PCRrm_bam.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input.sbam} 1>{log.out} 2>{log.err}"


rule get_ran_CG:
    input:
        refG=refG
    output:
        pozF=temp("aux_files/"+re.sub('.fa*','.poz.gz',os.path.basename(refG))),
        ranCG=os.path.join("aux_files",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
    log:
        err="aux_files/logs/get_ran_CG.err"
    threads: 1
    conda: mCtCondaEnvironment
    shell: 'set +o pipefail; ' + os.path.join(workflow_tools,'methylCtools') + " fapos {input.refG}  " + re.sub('.gz','',"{output.pozF}") + ';cat '+ re.sub('.gz','',"{output.pozF}") +' | grep "+" - | shuf | head -n 1000000 | awk \'{{print $1, $5, $5+1, $6, $8}}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + "{output.ranCG} 2>{log.err}"
        


rule calc_Mbias:
    input:
        refG=refG,
        rmDupBam="bams/{sample}.PCRrm.bam",
        sbami="bams/{sample}.PCRrm.bam.bai"
    output:
        mbiasTXT="QC_metrics/{sample}.Mbias.txt"
    log:
        out="QC_metrics/logs/{sample}.calc_Mbias.out"
    threads: nthreads
    conda: CondaEnvironment
    shell: "MethylDackel mbias {input.refG} {input.rmDupBam} {output.mbiasTXT} -@ {threads} 1>{log.out} 2>{output.mbiasTXT}"


if intList:
    rule depth_of_cov:
        input:
            refG=refG,
            rmDupBam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            ranCG=os.path.join("aux_files",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG))),
            intList=intList
        output:
            outFileList=calc_doc(intList,True)
        params:
            tempdir=tempdir,
            auxdir=os.path.join(wdir,"aux_files"),
            OUTlist=lambda wildcards,output: [w.replace('.sample_summary', '') for w in output.outFileList],
            OUTlist0=lambda wildcards,output: [w.replace('.sample_summary', '') for w in output.outFileList][0],
            OUTlist1=lambda wildcards,output: [w.replace('.sample_summary', '') for w in output.outFileList][1],
            auxshell=lambda wildcards,input,output: ';'.join("gatk -Xmx30g -Djava.io.tmpdir="+ tempdir +" -T DepthOfCoverage -R "+ input.refG +" -o "+ oi +" -I " + input.rmDupBam + " -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample -L " + bi for oi,bi in zip([w.replace('.sample_summary', '') for w in output.outFileList][2:],input.intList))
        log:
            err="QC_metrics/logs/{sample}.depth_of_cov.err",
            out="QC_metrics/logs/{sample}.depth_of_cov.out"
        threads: 1
        conda: CondaEnvironment
        shell: "gatk -Xmx30g -Djava.io.tmpdir={params.tempdir} -T DepthOfCoverage -R {input.refG} -o {params.OUTlist0} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample ; gatk -Xmx30g -Djava.io.tmpdir={params.tempdir}  -T DepthOfCoverage -R {input.refG} -o {params.OUTlist1} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample -L {input.ranCG}; {params.auxshell} 1>{log.out} 2>{log.err}" 


else:
    rule depth_of_cov:
        input:
            refG=refG,
            rmDupBam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            ranCG=os.path.join("aux_files",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
        output:
            outFileList=calc_doc(intList,True)
        params:
            tempdir=tempdir,
            auxdir=os.path.join(wdir,"aux_files"),
            OUTlist=lambda wildcards,output: [w.replace('.sample_summary', '') for w in output.outFileList]
        log:
            err="QC_metrics/logs/{sample}.depth_of_cov.err",
            out="QC_metrics/logs/{sample}.depth_of_cov.out"
        threads: 1
        conda: CondaEnvironment
        shell: "gatk -Xmx30g -Djava.io.tmpdir={params.tempdir} -T DepthOfCoverage -R {input.refG} -o {params.OUTlist0} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample ; gatk -Xmx30g -Djava.io.tmpdir={params.tempdir}  -T DepthOfCoverage -R {input.refG} -o {params.OUTlist1} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample -L {input.ranCG} 1>{log.out} 2>{log.err}"

if not trimReads is None:
    rule downsample_reads:
        input:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        output:
            R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        log:
            err="FASTQ_downsampled/logs/{sample}.downsample_reads.err",
            out="FASTQ_downsampled/logs/{sample}.downsample_reads.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: """
                seqtk sample -s 100 {input.R1cut} 5000000 | pigz -p {threads} -9 > {output.R1downsampled}  
                seqtk sample -s 100 {input.R2cut} 5000000 | pigz -p {threads} -9 > {output.R2downsampled}
                1>{log.out} 2>{log.err} 
               """

    rule conv_rate:
        input:
            R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        output:
            R12cr="QC_metrics/{sample}.conv.rate.txt"
        params:
            read_root=lambda wildcards: "FASTQ_downsampled/"+wildcards.sample
        log:
            err="QC_metrics/logs/{sample}.conv_rate.err",
            out="QC_metrics/logs/{sample}.conv_rate.out"
        threads: 1
        shell: os.path.join(workflow_tools,'conversionRate_KS.sh ')+ "{params.read_root} {output.R12cr} 1>{log.out} 2>{log.err}"

else:
    rule downsample_reads:
        input:
            R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2="FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        log:
            err="FASTQ_downsampled/logs/{sample}.downsample_reads.err",
            out="FASTQ_downsampled/logs/{sample}.downsample_reads.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: """
                seqtk sample -s 100 {input.R1} 5000000 | pigz -p {threads} -9 > {output.R1downsampled}  
                seqtk sample -s 100 {input.R2} 5000000 | pigz -p {threads} -9 > {output.R2downsampled}
                1>{log.out} 2>{log.err} 
               """
    
    rule conv_rate:
        input:
            R1="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        output:
            R12cr="QC_metrics/{sample}.conv.rate.txt"
        log:
            err="QC_metrics/logs/{sample}.conv_rate.err",
            out="QC_metrics/logs/{sample}.conv_rate.out"
        threads: 1
        shell: os.path.join(workflow_tools,'conversionRate_KS.sh ')+ "FASTQ_downsampled/{sample} {output.R12cr} 1>{log.out} 2>{log.err}"

rule get_flagstat:
    input:
        rmDupbam="bams/{sample}.PCRrm.bam"
    output:
        fstat="QC_metrics/{sample}.flagstat"
    log:
        err="QC_metrics/logs/{sample}.get_flagstat.err"        
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools flagstat {input.rmDupbam} > {output.fstat} 2>{log.err}" 

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
    log:
        err="QC_metrics/logs/produce_report.err",
        out="QC_metrics/logs/produce_report.out"
    conda: RmdCondaEnvironment 
    threads: 1
    shell: "cp -v " + os.path.join(workflow_tools,"WGBS_QC_report_template.Rmd")+ " " + os.path.join("aux_files", "WGBS_QC_report_template.Rmd") + ';Rscript -e "rmarkdown::render(\''+os.path.join(wdir,"aux_files", "WGBS_QC_report_template.Rmd")+'\', params=list(QCdir=\'"' + os.path.join(wdir,"QC_metrics") +'"\' ), output_file =\'"'+ os.path.join(wdir,"QC_metrics",'QC_report.pdf"\'')+')"' + " 1>{log.out} 2>{log.err}"


if mbias_ignore=="auto":
    rule methyl_extract:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=refG,
            mbiasTXT="QC_metrics/{sample}.Mbias.txt"     
        output:
            methTab="methXT/{sample}_CpG.bedGraph"
        params:
            OUTpfx=lambda wildcards,output: os.path.join(wdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="methXT/logs/{sample}.methyl_extract.err",
            out="methXT/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CondaEnvironment
        shell: "mi=$(cat {input.mbiasTXT} | sed 's/Suggested inclusion options: //' );MethylDackel extract  -o {params.OUTpfx} -q 10 -p 20 $mi --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(wdir,"{input.rmDupbam}") + " 1>{log.out} 2>{log.err}"
            

else:
    rule methyl_extract:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=refG
        output:
            methTab="methXT/{sample}_CpG.bedGraph"
        params:
            mbias_ignore=mbias_ignore,
            OUTpfx=lambda wildcards,output: os.path.join(wdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="methXT/logs/{sample}.methyl_extract.err",
            out="methXT/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CondaEnvironment
        shell: "MethylDackel extract  -o {params.OUTpfx} -q 10 -p 20 {params.mbias_ignore} --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(wdir,"{input.rmDupbam}") + " 1>{log.out} 2>{log.err}"

if blackList is None:
    rule CpG_filt:
        input:
            methTab="methXT/{sample}_CpG.bedGraph"
        output:
            tabFilt="methXT/{sample}.CpG.filt2.bed"
        params:
            methDir=os.path.join(wdir,"methXT"),
            OUTtemp=lambda wildcards,input: os.path.join(wdir,re.sub('_CpG.bedGraph','.CpG.filt.bed',input.methTab))
        log:
            err="methXT/logs/{sample}.CpG_filt.err",
            out="methXT/logs/{sample}.CpG_filt.out"
        threads: 1
        conda: CondaEnvironment
        shell: "Rscript --no-save --no-restore " + os.path.join(workflow_tools,'WGBSpipe.POM.filt.R ') + "{params.methDir} {input.methTab};mv -v {params.OUTtemp} {output.tabFilt} 1>{log.out} 2>{log.err}"

else:
    rule CpG_filt:
        input:
            methTab="methXT/{sample}_CpG.bedGraph",
            blackListF=blackList
        output:
            tabFilt="methXT/{sample}.CpG.filt2.bed"
        params:
            methDir=os.path.join(wdir,"methXT"),
            OUTtemp=lambda wildcards,input: os.path.join(wdir,re.sub('_CpG.bedGraph','.CpG.filt.bed',input.methTab))
        log:
            err="methXT/logs/{sample}.CpG_filt.err",
            out="methXT/logs/{sample}.CpG_filt.out"
        threads: 1
        conda: CondaEnvironment
        shell: "Rscript --no-save --no-restore " + os.path.join(workflow_tools,'WGBSpipe.POM.filt.R ') + "{params.methDir} {input.methTab};bedtools intersect -v -a {params.OUTtemp} -b {input.blackListF} > {output.tabFilt} 1>{log.out} 2>{log.err}"


if sampleInfo:
    rule CpG_stats:
        input: expand("methXT/{sample}.CpG.filt2.bed",sample=samples)
        output:
            RDatAll='singleCpG_stats_limma/singleCpG.RData',
            Limdat='singleCpG_stats_limma/limdat.LG.RData',
            MetIN='singleCpG_stats_limma/metilene.IN.txt'
        params:
            statdir=os.path.join(wdir,'singleCpG_stats_limma'),
            sampleInfo=sampleInfo
        log:
            err="singleCpG_stats_limma/logs/CpG_stats.err",
            out="singleCpG_stats_limma/logs/CpG_stats.out"
        threads: 1
        conda: CondaEnvironment
        shell: "Rscript --no-save --no-restore " + os.path.join(workflow_tools,'WGBSpipe.singleCpGstats.limma.R ') + "{params.statdir} {params.sampleInfo} "  + os.path.join(wdir,"methXT") + " 1>{log.out} 2>{log.err}"

    rule run_metilene:
        input:
            MetIN='singleCpG_stats_limma/metilene.IN.txt',
            sampleInfo=sampleInfo
        output:
            MetBed='metilene_out/singleCpG.metilene.bed'
        params:
            DMRout=os.path.join(wdir,'metilene_out')
        log:
            err="metilene_out/logs/run_metilene.err"
        threads: nthreads
        conda: CondaEnvironment
        shell: 'metilene -a ' + list(set(pandas.read_table(sampleInfo)['Group']))[0] + ' -b ' + list(set(pandas.read_table(sampleInfo)['Group']))[1] + " -t {threads} {input.MetIN} | sort -k 1,1 -k2,2n > {output.MetBed}" + " 2>{log.err}"


    rule get_CG_metilene:
        input:
            refG=refG,
            MetBed='metilene_out/singleCpG.metilene.bed',
            pozF="aux_files/"+re.sub('.fa*','.poz.gz',os.path.basename(refG))
        output:
            imdF=temp("aux_files/"+re.sub('.fa*','.CpG.bed',os.path.basename(refG))),
            MetCG=os.path.join("aux_files",re.sub('.fa','.metilene.CpGlist.bed',os.path.basename(refG)))
        params:
            auxdir=os.path.join(wdir,"aux_files")            
        log:
            err="aux_files/logs/get_CG_metilene.err"
        threads: 1
        conda: CondaEnvironment
        shell: 'grep "+"' + " {input.pozF} "+ ' | awk \'{{print $1, $5, $5+1, $6, $8}}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + "{output.imdF};"+ "bedtools intersect -wa -a {output.imdF} -b {input.MetBed} > {output.MetCG}  2>{log.err};sleep 300"
            

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
            DMRout=os.path.join(wdir,'metilene_out')
        log:
            err="metilene_out/logs/cleanup_metilene.err",
            out="metilene_out/logs/cleanup_metilene.out"
        threads: 1
        conda: CondaEnvironment
        shell: 'Rscript --no-save --no-restore ' + os.path.join(workflow_tools,'WGBSpipe.metilene_stats.limma.R ') + "{params.DMRout} " + os.path.join(wdir,"{input.MetBed}") +' ' + os.path.join(wdir,"{input.MetCG}") + ' ' + os.path.join(wdir,"{input.Limdat}") + " {input.sampleInfo} {input.refG} 1>{log.out} 2>{log.err}" 


if intList:
    rule get_CG_per_int:
        input:
            intList=intList,
            refG=refG,
            imdF="aux_files/"+re.sub('.fa*','.CpG.bed',os.path.basename(refG))
        output:
            outList=run_int_aggStats(intList,False) 
        log:
            err="aux_files/logs/get_CG_per_int.err"
        params:
            auxshell=lambda wildcards,input,output: ';'.join(["bedtools intersect -wa -a "+ input.imdF + " -b " + bli + ' > ' + oli  for bli,oli in zip(input.intList,output.outList) ])+';sleep 300'
        threads: 1
        conda: CondaEnvironment
        shell:"{params.auxshell} 2>{log.err}"


    if sampleInfo:
        rule intAgg_stats:
            input:
                Limdat='singleCpG_stats_limma/limdat.LG.RData',
                intList=intList,
                refG=refG,
                sampleInfo=sampleInfo,
                aux=run_int_aggStats(intList,False)
            output:
                outFiles=run_int_aggStats(intList,sampleInfo)
            params:
                auxshell=lambda wildcards,input:';'.join(['Rscript --no-save --no-restore ' + os.path.join(workflow_tools,'WGBSpipe.interval_stats.limma.R ') + os.path.join(wdir,'aggregate_stats_limma ') + li +' '+ aui +' ' + os.path.join(wdir,input.Limdat) + ' '  + input.sampleInfo  for li,aui in zip(intList,[os.path.join(wdir,"aux_files",re.sub('.fa',re.sub('.bed','.CpGlist.bed',os.path.basename(x)),os.path.basename(refG))) for x in intList])])
            log:
                err="aggregate_stats_limma/logs/intAgg_stats.err",
                out="aggregate_stats_limma/logs/intAgg_stats.out"
            threads: 1
            conda: CondaEnvironment
            shell: "{params.auxshell} 1>{log.out} 2>{log.err}"


 
