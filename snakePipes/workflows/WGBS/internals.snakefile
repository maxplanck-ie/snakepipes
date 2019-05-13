import glob
import os
import sys

## Main variables ##############################################################


### Functions ##################################################################


### Variable defaults ##########################################################
## trim
fastq_dir = "FASTQ"
if trim:
    fastq_indir_trim = "FASTQ"
    if trim_prg == "trimgalore":
        fastq_dir = "FASTQ_TrimGalore"
    elif trim_prg == "cutadapt":
        fastq_dir = "FASTQ_Cutadapt"

### Initialization #############################################################

# Disable trimming if BAM files are input
if fromBam:
    trim = False
    fastqc = False

##check if genes_bed is available
if not 'genes_bed' in globals():
    genes_bed='NA'

if not fromBam:
    infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*{}'.format(ext))))
    samples = cf.get_sample_names(infiles, ext, reads)
    paired = cf.is_paired(infiles, ext, reads)

    if not samples:
        sys.exit("Error! NO samples found in dir {}!!!\n".format(str(indir or '')))

    if not paired:
        sys.exit("Error! Paired-end samples not detected. "
                 "WGBS workflow currently works only with paired-end samples {}!!!\n".format(str(indir or '')))
else:
    infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*{}'.format(bam_ext))))
    samples = cf.get_sample_names_bam(infiles, bam_ext)

if not samples:
    sys.exit("\n  Error! NO samples found in dir {}!!!\n".format(str(indir or '')))


def getGroups(sampleSheet):
    """
    Given a sample sheet, return a tuple of (group1, group2) for use with metilene

    Mut and Treatment come last, otherwise the lexographic order is returned
    """
    groups = set()
    conditionIdx = None
    f = open(sampleSheet)
    line = f.readline()
    cols = line.strip().split()
    if "condition" not in cols:
        conditionIdx = 2
    else:
        conditionIdx = cols.index("condition")
    for line in f:
        cols = line.strip().split()
        groups.add(cols[conditionIdx])
    f.close()

    if "Mut" in groups:
        return (list(groups - set(["Mut"]))[0], "Mut")
    if "Treatment" in groups:
        return (list(groups - set(["Treatment"]))[0], "Treatment")
    return sorted(list(groups))[:2]
