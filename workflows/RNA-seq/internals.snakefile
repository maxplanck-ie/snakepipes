import glob
import os
import subprocess

#print(vars(workflow))
#print("CONFIG:", config)
# print("overwrite_configfile:", workflow.overwrite_configfile)

## Load defaults.yaml with default parameters
with open(os.path.join(workflow.basedir, "defaults.yaml"), "r") as f:
    defaults = yaml.load(f)


##print("\n--- defaults.yaml --------------------------------------------------------------")
##for k,v in sorted(defaults.items()):
##    print("{}: {}".format(k,v))
##print


## Load .config.yaml with basic configuration settings from wrapper ############
def merge_dicts(x, y):
    z = x.copy()
    z.update(y)
    return(z)

try:
    if workflow.overwrite_configfile:
        pass # automatically uses snakemake config file (--configfile)
    else:
        ## or use config file from the output dir
        with open(os.path.join(workflow.workdir_init, "config.yaml"), "r") as f:
            config = yaml.load(f)
    config = merge_dicts(defaults, config)
except:
    config = defaults


## Main variables ##############################################################
maindir = os.path.dirname(os.path.dirname(workflow.basedir))
verbose = config["verbose"]

if verbose:
    print("\n--- config ----------------------------------------------------------------")
    for k,v in sorted(config.items()):
        print("{}: {}".format(k,v))
    print()


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

indir = config["indir"]
outdir = config["outdir"]
reads = config["reads"]
library_type = config["library_type"]
ext = config["ext"]
downsample = config["downsample"]
genome = config["genome"]
fragment_length = config["fragment_length"]
filter_annotation = config["filter_annotation"]
if not filter_annotation:
    filter_annotation = "''"
salmon_index_options = config["salmon_index_options"]
trim_galore_opts = config["trim_galore_opts"]
fastqc = config["fastqc"]

if config["trim"] == "trimgalore":
    trim = True
    fastq_dir = "FASTQ_TrimGalore"
    fastq_indir_trim = "FASTQ"
elif config["trim"] == "cutadapt":
    trim = True
    fastq_dir = "FASTQ_Cutadapt"
    fastq_indir_trim = "FASTQ"
else:
    trim = False
    fastq_dir = "FASTQ"


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
