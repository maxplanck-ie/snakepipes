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
strains = list(map(str.strip, re.split(",|;", config["strains"])))


infiles = sorted(glob.glob(os.path.join(str(indir or ""), "*" + ext)))
samples = cf.get_sample_names(infiles, ext, reads)
pairedEnd = cf.is_paired(infiles, ext, reads)
del infiles


# reference genomes: diploid_genome, nmasked, reference
# currently only the case for diploid_genome is covered
REFERENCES = ["diploid_genome"]

# possible phasetypes of contacts
PHASETYPES = [strains[0], strains[1], "unphased", "trans"]

PHASEFILTER = [
    '(phase1=="0") and (phase2=="0")',
    '(phase1=="1") and (phase2=="1")',
    '(phase1==".") or (phase2==".")',
    '(phase1!=phase2) and (phase1!=".") and (phase2!=".") and (phase1!="!") and (phase2!="!")',
]

PHASEDIC = dict(map(lambda i, j: (i, j), PHASETYPES, PHASEFILTER))
