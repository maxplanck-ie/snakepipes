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
    


def get_Read_Group(INfile1):
    with io.TextIOWrapper(gzip.open(INfile1, 'r')) as f:
        file_content = f.readline().strip()
    PL=re.sub('@','',file_content).split(":")[0]
    PU=re.sub('@','',file_content).split(":")[2]
    RG='@RG"\t"ID:1"\t"SM:{sample}"\t"LB:{sample}"\t"PL:'+PL+'"\t"PU:'+PU
    
    return RG


def get_mbias_auto():
    with open(os.path.join("QC_metrics","{sample}.Mbias.txt", 'r') as f:
            first_line = f.readline()
        m_ignore=re.sub('Suggested inclusion options: ','',first_line).strip('\n')
    return m_ignore



