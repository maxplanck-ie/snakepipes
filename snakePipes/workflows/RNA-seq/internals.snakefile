import glob
import os
import subprocess
import re


## Main variables ##############################################################


### Functions ##################################################################


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
if not fromBam:
    infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
    samples = cf.get_sample_names(infiles,ext,reads)

    paired = cf.is_paired(infiles,ext,reads)

    if not paired:
        reads = [""]
else:
    infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*'+bam_ext)))
    samples = cf.get_sample_names_bam(infiles,bam_ext)


if sampleSheet:
    cf.check_sample_info_header(sampleSheet)

if sampleSheet and not cf.check_replicates(sampleSheet):
    print("\nWarning! Sleuth cannot be invoked without replicates! Only DESeq2 is used...\n")
