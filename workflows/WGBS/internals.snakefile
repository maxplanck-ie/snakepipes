import os
import re
import glob
import subprocess
import zipfile
import pandas
import gzip
import io
import tempfile
import shutil

def calc_cutThd (zipL,fqin,logobject,wdir):
    rNcutL=[]
    with open(logobject,"w") as lo:
        for zipi in zipL: 
            zf=os.path.basename(zipi)
            print("Processing zipped fastqc file:" + zf,file=lo)
            if not os.path.exists(os.path.join(wdir,"FastQC_In",re.sub('\.zip','',zf))):
                with zipfile.ZipFile(zipi, "r") as z:
                        z.extractall(path=os.path.join(wdir,"FastQC_In"))
            fqtxt=os.path.join(wdir,"FastQC_In",re.sub('\.zip','',zf),'fastqc_data.txt')
            print('Currently processing :'+ fqtxt,file=lo)
            os.chdir(os.path.join(wdir,"FastQC_In",re.sub('\.zip','',zf)))
            subprocess.check_output(['csplit', '-z' , fqtxt , '/>>/','{*}'])
            with open(fqtxt,'r') as file:
                line=file.readline().strip()
            if '0.11.2' in line or '0.11.6' in line:
                NTconTab=pandas.read_table(os.path.join(os.getcwd(), 'xx09'), sep='\t',skiprows=1,header=0,names=['Index','G','A','T','C'],dtype={'Index':'object','G':'float64','A':'float64','T':'float64','C':'float64'},engine='c')
            else:
                NTconTab=pandas.read_table(os.path.join(os.getcwd(), 'xx09'), sep='\t',skiprows=1,header=0,names=['Index','G','A','T','C'],dtype={'Index':'object','G':'float64','A':'float64','T':'float64','C':'float64'},engine='c')
                print('Check fastqc version',file=lo)
            difftab=NTconTab.set_index('Index').diff(periods=-1)
            difftabA=difftab.abs()
            maxv=difftabA.idxmax(axis=0)
            maxv=maxv.values.astype(int)
            rNmax=list(difftabA.index)
            rNcut=rNmax[(maxv.max()-1)]
            rNcutL.append(str(rNcut)) ##
            print(NTconTab.head(n=10),file=lo)
            print(difftab.head(n=10),file=lo)
            print('Maximal absolute difference per nucleotide :',file=lo)
            print(difftabA.max(axis=0),file=lo)
            print('Index of diffmax :',file=lo)
            print(difftabA.idxmax(axis=0),file=lo)
            print('Index of the maximal difference :',file=lo)
            print(maxv.max(),file=lo)
            print('Number of nucleotides for 5prime trimming :' + rNcut,file=lo)
            os.getcwd()
    zipLre=[ re.sub('_fastqc.zip','.fastq.gz',x ) for x in zipL ]    
    cutThdRes=OrderedDict(zip(zipLre, rNcutL))
    ctr1=filter(lambda x:'_R1.fastq.gz' in x, cutThdRes.keys())
    ctr2=filter(lambda x:'_R2.fastq.gz' in x, cutThdRes.keys())
    cutThdRes_R1=[ cutThdRes[x] for x in ctr1 ]
    cutThdRes_R2=[ cutThdRes[x] for x in ctr2 ]
    cutThdL=zip(cutThdRes_R1,cutThdRes_R2)
    return cutThdL
    


def cut_reads_auto(inFile):
    #prepare threshold values
    with open(inFile,'r') as ctfile:
        for line in ctfile:
            fields = line.strip().split()
            ct1=fields[0]
            ct2=fields[1]
    bshcmd=cutadapt_module_path +" cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5 -j {threads}" + ' -u ' + ct1 + ' -U ' + ct2 + " -o {output.R1cut} -p {output.R2cut} {input.R1} {input.R2}" 
    return bshcmd

def cut_reads_user(trimThreshold,trimOtherArgs,nextera):
    adapterSeq = "AGATCGGAAGAGC"
    if nextera:
        adapterSeq = "CTGTCTCTTATA"
    bshcmd = "{} cutadapt -a {} -A {} -q {} -m 30 -j {threads} {} -o {} -p {} {} {} ".format(cutadapt_module_path,
                                                                            adapterSeq,
                                                                            adapterSeq,
                                                                            trimThreshold,
                                                                            trimOtherArgs,
                                                                            {output.R1cut},
                                                                            {output.R2cut},
                                                                            {input.R1},
                                                                            {input.R2})
    return(bshcmd)



def bMeth_map_reads(INfile1,INfile2,crefG):
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(INfile1))
    with io.TextIOWrapper(gzip.open(INfile1, 'r')) as f:
        file_content = f.readline().strip()
    PL=re.sub('@','',file_content).split(":")[0]
    PU=re.sub('@','',file_content).split(":")[2]
    RG='@RG"\t"ID:1"\t"SM:'+read_root+'"\t"LB:'+read_root+'"\t"PL:'+PL+'"\t"PU:'+PU
    # There no point in using more than a few sort threads, it just uses excess memory
    sortThreads = nthreads
    if nthreads > 4:
        sortThreads = 4
    mapcmd=bwameth_path + " bwameth.py --threads  {threads}  --read-group "+ RG + ' --reference ' + crefG + ' ' + INfile1 + ' ' + INfile2 + ' | ' + os.path.join(samtools_path,'samtools') + ' sort -T ' + tempfile.mkdtemp(suffix='',prefix=read_root,dir='/data/extended') + ' -m 3G -@ ' + str(sortThreads)  + " -o {output.sbam}"
    return mapcmd



def mCT_get_ranCG():
    #pozF=os.path.join("aux_files",re.sub('.fa*','.poz.gz',os.path.basename(input.refG)))
    #'source activate NGSpy2.7; ' +' ;source deactivate'
    cmd_from0='module load methylCtools;set +o pipefail; ' + os.path.join(mCT_path,'methylCtools') + ' fapos ' + "{input.refG}" + ' ' + re.sub('.gz','',"{output.pozF}") + ';cat '+ re.sub('.gz','',"{output.pozF}") +' | grep "+" - | shuf | head -n 1000000 | awk \'{{print $1, $5, $5+1, $6, $8}}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + "{output.ranCG}" 
    cmd_from_poz='zcat ' +  "{output.pozF}" + ' | grep "+" - | shuf | head -n 1000000 | awk \'{{print $1, $5, $5+1, $6, $8}}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + "{output.ranCG}"
    if os.path.exists("{output.pozF}"):
        shcmd=cmd_from_poz
    else:
        shcmd=cmd_from0
    return shcmd



def BS_doc_XT(outListAuto,intList,refG):
    OUTlist=[w.replace('.sample_summary', '') for w in outListAuto]
    OUTlist2=OUTlist[2:]
    WG_mean_cmd='java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATK_path,'GenomeAnalysisTK.jar')+" -R {input.refG} -T DepthOfCoverage -o " + str(OUTlist[0]) + " -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample "
    CG_cmd='java -Xmx50g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATK_path,'GenomeAnalysisTK.jar')+" -R {input.refG} -T DepthOfCoverage -o " + str(OUTlist[1]) + " -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -omitIntervals -mmq 10 --partitionType sample -L " + os.path.join("{params.auxdir}",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
    cmd_Xi=['java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATK_path,'GenomeAnalysisTK.jar')+" -R {input.refG} -T DepthOfCoverage -o " + oi + " -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -omitIntervals -mmq 10 --partitionType sample -L " + bi for oi,bi in zip(OUTlist2,intList)]
    cmd_all=cmd_Xi
    cmd_all[0:0]=[CG_cmd]
    cmd_all[0:0]=[WG_mean_cmd]
    cmd_all_str=';'.join(cmd_all)
    return cmd_all_str



def BS_doc(outListAuto,refG):
    OUTlist=[w.replace('.sample_summary', '') for w in outListAuto]
    WG_mean_cmd='java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATK_path,'GenomeAnalysisTK.jar')+" -R {input.refG} -T DepthOfCoverage -o " + str(OUTlist[0]) + " -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample "
    CG_cmd='java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATK_path,'GenomeAnalysisTK.jar')+" -R {input.refG} -T DepthOfCoverage -o " + str(OUTlist[1]) + " -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -omitIntervals -mmq 10 --partitionType sample -L " + os.path.join("{params.auxdir}",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
    cmd_all=[WG_mean_cmd,CG_cmd]
    cmd_all_str=';'.join(cmd_all)
    return cmd_all_str


def BS_QC_rep():
    shutil.copyfile(os.path.join(workflow_tools,"WGBS_QC_report_template.Rmd"), os.path.join("aux_files", "WGBS_QC_report_template.Rmd"))
    cmd = os.path.join(R_path, 'Rscript -e "rmarkdown::render(\''+os.path.join(wdir,"aux_files", "WGBS_QC_report_template.Rmd")+'\', params=list(QCdir=\'"' + os.path.join(wdir,"QC_metrics") +'"\' ), output_file =\'"'+ os.path.join(wdir,"QC_metrics",'QC_report.pdf"\'')+')"')
    return cmd


def methXT_POM(inFile,mbias_action,wdir,OutFile):
    read_root=re.sub('.PCRrm.bam','',os.path.basename(inFile)) 
    OUTpfx=os.path.join(wdir,re.sub('_CpG.bedGraph','',OutFile))
    if len(mbias_action) < 3:
        m_ignore=','.join([mbias_ignore]*4)
        POM_cmd=os.path.join(POM_path,'MethylDackel') + ' extract  -o ' + OUTpfx + ' -q 10 -p 20 --nOT ' + m_ignore + ' --nOB  ' + m_ignore + " --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(wdir,"{input.rmDupbam}") 
    elif mbias_action=="auto":
        with open(os.path.join("QC_metrics",read_root)+'.Mbias.txt', 'r') as f:
            first_line = f.readline()
        m_ignore=re.sub('Suggested inclusion options: ','',first_line).strip('\n')
        POM_cmd=os.path.join(POM_path,'MethylDackel') + ' extract  -o ' + OUTpfx + ' -q 10 -p 20 ' + m_ignore + " --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(wdir,inFile) 
    elif len(mbias_action) > 4:
        POM_cmd=os.path.join(POM_path,'MethylDackel') + ' extract  -o ' + OUTpfx + ' -q 10 -p 20 ' + mbias_action + " --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(wdir,inFile) 
    return POM_cmd



def filt_POM(inFile,blackList):
    read_root=re.sub('_CpG.bedGraph','',os.path.basename(inFile))
    Rfilt_cmd=os.path.join(R_path,'Rscript') +' --no-save --no-restore ' + os.path.join(workflow_tools,'WGBSpipe.POM.filt.R ') + "{params.methDir} {input.methTab} {params.Rlib}" 
    if blackList is None:
        mv_cmd='mv -v '+ re.sub('_CpG.bedGraph','.CpG.filt.bed',inFile) + ' ' + re.sub('_CpG.bedGraph','.CpG.filt2.bed',inFile) 
        cmd_all=';'.join([Rfilt_cmd,mv_cmd])
    else:
        SNP_filt=bedtools_module_path +' bedtools' + ' intersect -v -a' + re.sub('_CpG.bedGraph','.CpG.filt.bed',inFile) + ' -b ' + blackList + ' > ' + re.sub('_CpG.bedGraph','.CpG.filt2.bed',inFile)
        cmd_all=';'.join([Rfilt_cmd,SNP_filt])
    return cmd_all



def DMR_metilene(sampleInfo):
    met_cmd=os.path.join(metilene_path,'metilene')+ ' -a ' + list(set(pandas.read_table(sampleInfo)['Group']))[0] + ' -b ' + list(set(pandas.read_table(sampleInfo)['Group']))[1] + " -t {threads} {input.MetIN} | sort -k 1,1 -k2,2n > {output.MetBed}" 
    return met_cmd



def mCT_get_CpGxInt(pozF,outList,bedList,wdir):
    #pozF=os.path.join(wdir,"aux_files",re.sub('.fa*','.poz',os.path.basename(refG)))
    imdF=os.path.join(wdir,"aux_files",re.sub('.poz','.CpG.bed',os.path.basename(pozF)))
    if not os.path.exists(imdF):
        cmd0='grep "+"' + " {input.pozF} "+ ' | awk \'{{print $1, $5, $5+1, $6, $8}}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + imdF
        cmd1=[bedtools_module_path +' bedtools intersect -wa -a ' + imdF + ' -b ' + bli + ' > ' + oli for bli,oli in zip(bedList,outList) ]
        cmd_all=cmd1
        cmd_all[0:0]=[cmd0]
        cmd_all_str=';'.join(cmd_all)+';sleep 300'
    else:
        cmd1=[bedtools_module_path +' bedtools intersect -wa -a ' + imdF + ' -b ' + bli + ' > ' + oli for bli,oli in zip(bedList,outList) ]
        cmd_all=cmd1
        cmd_all_str=';'.join(cmd_all)+';sleep 300'
    return cmd_all_str



def clean_up_metilene():
    cmd=os.path.join(R_path,'Rscript') +' --no-save --no-restore ' + os.path.join(workflow_tools,'WGBSpipe.metilene_stats.limma.R ') + "{params.DMRout} " + os.path.join(wdir,"{input.MetBed}") +' ' + os.path.join(wdir,"{input.MetCG}") + ' ' + os.path.join(wdir,"{input.Limdat}") + " {input.sampleInfo} {input.refG}" + ' "' + bedtools_module_path+'" ' + "{params.Rlib}" 
    return cmd



def int_stats_limma(intList):
    auxList=[os.path.join("{params.auxdir}",re.sub('.fa',re.sub('.bed','.CpGlist.bed',os.path.basename(x)),os.path.basename(refG))) for x in intList]
    cmd_all=[os.path.join(R_path,'Rscript') +' --no-save --no-restore ' + os.path.join(workflow_tools,'WGBSpipe.interval_stats.limma.R ') + "{params.aggStatdir} " + li +' '+ aui +' ' + os.path.join(wdir,"{input.Limdat}") + " {input.sampleInfo} {params.Rlib}" for li,aui in zip(intList,auxList)]
    cmd_all_str=';'.join(cmd_all)
    return cmd_all_str


