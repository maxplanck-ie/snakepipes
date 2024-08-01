
import glob
import os

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


## genome names for allele-sp mapping
strains = list(map(str.strip, re.split(',|;', config["strains"])))


infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*'+ext)))
samples = cf.get_sample_names(infiles,ext,reads)
pairedEnd = cf.is_paired(infiles,ext,reads)
del infiles
