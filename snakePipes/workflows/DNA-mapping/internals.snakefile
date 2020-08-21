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
if infiles == []:
    sys.exit("Error! Samples extnesion in {} are not {}. "
             "Please change the extensions to it or update the config.yaml file "
             "with your desired extension.".format(indir,ext))
samples = cf.get_sample_names(infiles,ext,reads)

pairedEnd = cf.is_paired(infiles,ext,reads)

del infiles

if not samples:
    sys.exit("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")

fromBAM = None

idxRange = 1
if pairedEnd:
    idxRange = 2

# clean up filtered_bam if needed appropriate
os.makedirs(outdir, exist_ok=True)
filt = ""
if dedup:
    filt += "-F 1024 "
    assert UMIDedup is False, "\nPlease use either --UMIDedup (UMI-tools dedup) or --dedup (via sambamba and samtools)!\n"
    "should be called!"
if properPairs:
    filt += "-f 2 "
if mapq is not None and mapq > 0:
    filt += "-q {} ".format(mapq)
filter_rules = os.path.join(outdir, "filter_rules")
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
