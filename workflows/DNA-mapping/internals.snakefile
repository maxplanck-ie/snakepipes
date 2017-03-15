import glob
import os
import subprocess
import re

#print(vars(workflow))
#print("CONFIG:", config)
# print("overwrite_configfile:", workflow.overwrite_configfile)


## Load defaults.yaml with default parameters ##################################

with open(os.path.join(workflow.basedir, "defaults.yaml"), "r") as f:
    defaults = yaml.load(f)

##print("\n--- defaults.yaml --------------------------------------------------------------")
##for k,v in sorted(defaults.items()):
##    print("{}: {}".format(k,v))
##print


## Load .config.yaml with basic configuration settings from wrapper
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


# When modifying the function update_filter(), double-check wether the rule
# samtools_filter has to be modified concordantly
def update_filter(samples):
    """
    Ensure that only the specified filters are applied
    If filtered BAM and sample.filter files exist already, check that they
    fully agree with the specified filtering options
    """
    # string with samtools view parameters for filtering
    filter = ""
    if dedup:
        filter += "-F 1024 "
    if properpairs:
        filter += "-f 2 "
    if mapq > 0:
        filter += "-q "+str(mapq)+" "

    for sample in samples:
        filtered_bam = os.path.join(outdir, "filtered_bam/"+sample+".filtered.bam")
        filtered_bai = filtered_bam+".bai"
        # filtered BAM file index sample.filtered.bam.bai exists already
        if os.path.isfile(filtered_bai):
            filter_file = filtered_bam.replace(".filtered.bam",".filter")
            # file sample.filter does not exist, i.e. sample.filtered.bam is
            # symlink to Bowtie2 output BAM file (without any filtering)
            if not os.path.isfile(filter_file):
                # remove symlink if filtering options have been specified
                if filter:
                    cmd = "rm -f "+filtered_bam+" "+filtered_bai
                    print("\nWARNING: Filtering options changed.\n"
                          "Removing files:", filtered_bam, filtered_bai)
                    subprocess.call(cmd, shell=True)
            # file sample.filter exists, i.e. filtering was done previously
            else:
                with open(filter_file, "r") as f:
                    old_filter = f.readline().strip()
                    current_filter = ("samtools view arguments: "+filter).strip()
                    # if specified filtering options differ from previous ones,
                    # then remove filtered files
                    if not current_filter == old_filter:
                        cmd = "rm -f "+filtered_bam+" "+filtered_bai+" "+filter_file
                        print("\nWARNING: Filtering options changed.\n"
                              "Removing files:", filtered_bam, filtered_bai, filter_file)
                        subprocess.call(cmd, shell=True)

    return


### Variable defaults ##########################################################

## Import from config into global name space! DANGEROUS!!!
for k,v in config.items():
    globals()[k] = v


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

# ensure that only the specified filters are applied to all files
# delete already filtered BAM files if they were generated with different
# filtering parameters in previous runs
update_filter(samples)
