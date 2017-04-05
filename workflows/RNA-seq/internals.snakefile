import glob
import os
import subprocess
import re

#print(vars(workflow))


## Main variables ##############################################################
maindir = os.path.dirname(os.path.dirname(workflow.basedir))
verbose = config["verbose"]


### Functions ##################################################################

def get_sample_names(infiles):
    """
    Get sample names without file extensions
    """
    s = []
    for x in infiles:
        x = os.path.basename(x).replace(ext,"")
        try:
            x = x.replace(reads[0],"").replace(reads[1],"")
        except:
            pass
        s.append(x)
    return(sorted(list(set(s))))


def is_paired(infiles):
    """
    Check for paired-end input files
    """
    paired = False
    infiles_dic = {}
    for infile in infiles:
        fname = os.path.basename(infile).replace(ext, "")
        m = re.match("^(.+)("+reads[0]+"|"+reads[1]+")$", fname)
        if m:
            ##print(m.group())
            bname = m.group(1)
            ##print(bname)
            if bname not in infiles_dic:
                infiles_dic[bname] = [infile]
            else:
                infiles_dic[bname].append(infile)
    if infiles_dic and max([len(x) for x in infiles_dic.values()]) == 2:
        paired = True
    # TODO: raise exception if single-end and paired-end files are mixed
    return(paired)


def convert_library_type (paired, from_library_type, from_prg, to_prg,
                            rscript=os.path.join(maindir, "shared", "tools", "library_type.R"),
                            tsv=os.path.join(maindir, "shared", "tools", "library_type.tsv")):
    """ Converts the library to e.g. from 2 (featureCounts) to RF (HISAT2) """
    if paired:
        lib_str = "PE"
    else:
        lib_str = "SE"

    cmd = ("Rscript {} {} {} {} {} {}".format(rscript, tsv, lib_str, from_library_type, from_prg, to_prg) )
    ##print("\n"+cmd)

    return( subprocess.check_output(cmd, shell=True).decode() )


## Variable defaults ##########################################################

print("\n--- config ---------------------------------------------------------------------")
for k,v in sorted(config.items()):
    globals()[k] = v    ## Import from config into global name space! DANGEROUS!!!
    if verbose:
        print("{}: {}".format(k,v))
print()


mode = list(map( str.strip, re.split(',|;', config["mode"]) ))

if not filter_annotation:
    filter_annotation = "''"

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
samples = get_sample_names(infiles)

paired = is_paired(infiles)

if not paired:
    reads = [""]

## rna-strandness for HISAT2
rna_strandness = convert_library_type(paired, library_type, "featureCounts", "HISAT2")
if rna_strandness == "NA":
    rna_strandness = ""
else:
    rna_strandness = "--rna-strandness "+rna_strandness

salmon_libtype = convert_library_type(paired, library_type, "featureCounts", "Salmon")
