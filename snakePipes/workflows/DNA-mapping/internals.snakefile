import glob
import os
import subprocess

## Main variables ##############################################################


### Functions ##################################################################


### Variable defaults ##########################################################
mode = list(map( str.strip, re.split(',|;', config["mode"]) ))
mode = [element.lower() for element in mode]
## genome names for allele-sp mapping
strains = list(map( str.strip, re.split(',|;', config["strains"]) ))
## trim
fastq_dir = "FASTQ"
if trim:
    fastq_indir_trim = "FASTQ"
    if trimmer == "trimgalore":
        fastq_dir = "FASTQ_TrimGalore"
    elif trimmer == "cutadapt":
        fastq_dir = "FASTQ_Cutadapt"
    elif trimmer == "fastp":
        fastq_dir = "FASTQ_fastp"


### Initialization #############################################################

infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*'+ext)))

pairedEnd = cf.is_paired(infiles,ext,reads)

samples = cf.get_sample_names(infiles,ext,reads)

del infiles

if not samples:
    print("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")
    exit(1)

fromBAM = None

if not fromBAM:
    cf.check_gz_reads(indir)

idxRange = 1
if pairedEnd:
    idxRange = 2
