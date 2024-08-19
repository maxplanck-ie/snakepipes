import glob
import os
import subprocess
import re
import sys


## Main variables ##############################################################


### Functions ##################################################################


## Variable defaults ##########################################################

mode = list(map(str.strip, re.split(',|;', config["mode"])))

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
if not fromBAM:
    infiles = sorted(glob.glob(os.path.join(indir, '*' + ext)))
    samples = cf.get_sample_names(infiles, ext, reads)
    pairedEnd = cf.is_paired(infiles, ext, reads)
else:
    infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*' + bamExt)))
    samples = cf.get_sample_names_bam(infiles, bamExt)

if sampleSheet:
    cf.check_sample_info_header(sampleSheet)
    isMultipleComparison = cf.isMultipleComparison(sampleSheet)

if sampleSheet and not cf.check_replicates(sampleSheet):
    sys.stderr.write("\nWarning! Sleuth cannot be invoked without replicates! Only DESeq2 is used...\n")

if not samples:
    print("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")
    exit(1)
