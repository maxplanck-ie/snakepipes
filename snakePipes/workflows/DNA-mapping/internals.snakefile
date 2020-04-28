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

idxRange = 1
if pairedEnd:
    idxRange = 2

# clean up filtered_bam if needed appropriate
os.makedirs(args.outdir, exist_ok=True)
filt = ""
if dedup:
    filt += "-F 1024 "
    assert args.UMIDedup is False, "\nPlease use either --UMIDedup (UMI-tools dedup) or --dedup (via sambamba and samtools)!\n"
    "should be called!"
if args.properPairs:
    filt += "-f 2 "
if mapq is not None and mapq > 0:
    filt += "-q {} ".format(args.mapq)
filter_rules = os.path.join(args.outdir, "filter_rules")
if os.path.exists(filter_rules):
    f = open(filter_rules)
    cont = f.read()
    f.close()
    if cont != filt:
        f = open(filter_rules, "w")
        f.write(filt)
        f.close()
else:
    f = open(filter_rules, "w")
    f.write(filt)
    f.close()
