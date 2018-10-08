import glob
import os
import subprocess
import re


## Main variables ##############################################################


### Functions ##################################################################

def check_replicates(sample_info_file):
    """
    return True if each condition has at least 2 replicates
    this check is eg. necessary for sleuth  
    """
    ret = subprocess.check_output(
            "cat "+sample_info_file+"| awk '/^\S*$/{next;}{if (NR==1){ col=0; for (i=1;i<=NF;i++) if ($i~\"condition\") col=i}; if (NR>1) print $col}' | sort | uniq -c | awk '{if ($1>1) ok++}END{if (NR>1 && ok>=NR) print \"REPLICATES_OK\"}'",
            shell=True).decode()

    if ret.find("REPLICATES_OK") >=0:
        return True
    else:
        return False


## Variable defaults ##########################################################

mode = list(map(str.strip, re.split(',|;', config["mode"])))

## genome names for allele-sp mapping
strains = list(map(str.strip, re.split(',|;', config["strains"])))

## trim
fastq_dir = "FASTQ"
if trim:
    fastq_indir_trim = "FASTQ"
    if trim_prg == "trimgalore":
        fastq_dir = "FASTQ_TrimGalore"
    elif trim_prg == "cutadapt":
        fastq_dir = "FASTQ_Cutadapt"


### Initialization #############################################################

infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = cf.get_sample_names(infiles,ext,reads)

paired = cf.is_paired(infiles,ext,reads)

if not paired:
    reads = [""]

if sample_info and not os.path.isfile(sample_info):
    print("ERROR: Cannot find sample info file! ("+sample_info+")\n")
    exit(1)

if sample_info and not cf.check_sample_info_header(sample_info):
    print("ERROR: Please use 'name' and 'condition' as column headers in sample info file! ("+sample_info+")\n")
    exit(1)

if sample_info and not check_replicates(sample_info):
    print("\nWarning! Sleuth cannot be invoked without replicates! Only DESeq2 is used...\n")
