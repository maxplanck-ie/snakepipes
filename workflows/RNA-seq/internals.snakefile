import glob
import os
import subprocess
import re


## Main variables ##############################################################

### Functions ##################################################################

def convert_library_type (paired, from_library_type, from_prg, to_prg,
                            rscript=os.path.join(maindir, "shared", "tools", "library_type.R"),
                            tsv=os.path.join(maindir, "shared", "tools", "library_type.tsv")):
    """ Converts the library to e.g. from 2 (featureCounts) to RF (HISAT2) """
    if paired:
        lib_str = "PE"
    else:
        lib_str = "SE"

    cmd = ("{}Rscript {} {} {} {} {} {}".format(R_path, rscript, tsv, lib_str, from_library_type, from_prg, to_prg) )
    ##print("\n"+cmd)

    return( subprocess.check_output(cmd, shell=True).decode() )


## Variable defaults ##########################################################

mode = list(map( str.strip, re.split(',|;', config["mode"]) ))
mode = [element.lower() for element in mode]

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

## rna-strandness for HISAT2
rna_strandness = convert_library_type(paired, library_type, "featureCounts", "HISAT2")
if rna_strandness == "NA":
    rna_strandness = ""
else:
    rna_strandness = "--rna-strandness "+rna_strandness

salmon_libtype = convert_library_type(paired, library_type, "featureCounts", "Salmon")
