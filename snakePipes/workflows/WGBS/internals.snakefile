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
    if trimmer == "trimgalore":
        fastq_dir = "FASTQ_TrimGalore"
    elif trimmer == "cutadapt":
        fastq_dir = "FASTQ_Cutadapt"
    elif trimmer == "fastp":
        fastq_dir = "FASTQ_fastp"

### Initialization #############################################################

# Disable trimming if BAM files are input
if fromBAM:
    trim = False
    aligner = "EXTERNAL_BAM"

##check if genes_bed is available
if not 'genes_bed' in globals():
    genes_bed='NA'

if not fromBAM:
    infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*{}'.format(ext))))
    samples = cf.get_sample_names(infiles, ext, reads)
    pairedEnd = cf.is_paired(infiles, ext, reads)

    if not samples:
        sys.exit("Error! NO samples found in dir {}!!!\n".format(str(indir or '')))

    if not pairedEnd:
        sys.exit("Error! Paired-end samples not detected. "
                 "WGBS workflow currently works only with paired-end samples {}!!!\n".format(str(indir or '')))
else:
    infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*{}'.format(bamExt))))
    samples = cf.get_sample_names_bam(infiles, bamExt)

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
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            if "condition" not in cols or "name" not in cols:
                sys.exit("ERROR: Please use 'name' and 'condition' as column headers in the sample info file ({})!\n".format(sample_info_file))
            conditionIdx = cols.index("condition")
            nCols = len(cols)
            continue
        elif idx == 1:
            # Sometimes there's a column of row names, which lack a header
            if len(cols) != nCols and len(cols) - 1 != nCols:
                sys.exit("ERROR: there's a mismatch between the number of columns in the header and body of {}!\n".format(sample_info_file))
            if len(cols) - 1 == nCols:
                conditionIdx += 1
        if idx > 0:        
            groups.add(cols[conditionIdx])
    f.close()

    if "Mut" in groups:
        return (list(groups - set(["Mut"]))[0], "Mut")
    if "Treatment" in groups:
        return (list(groups - set(["Treatment"]))[0], "Treatment")
    return sorted(list(groups))[:2]
