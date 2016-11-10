import glob
import os
import subprocess


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


# When modifying the function update_filter(), double-check if the rule
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

try:
    indir = config["indir"]
except:
    indir = os.getcwd()

try:
    outdir = config["outdir"]
except:
    outdir = os.getcwd()

## paired-end read name extension
try:
    reads = config["reads"]
except:
    reads = ["_R1", "_R2"]

## FASTQ file extension
try:
    ext = config["ext"]
except:
    ext = ".fastq.gz"

## downsampling - number of reads
try:
    downsample = int(config["downsample"])
except:
    downsample = None

try:
    genome = config["genome"]
except:
    genome = None

try:
    mate_orientation = config["mate_orientation"]
except:
    mate_orientation = "--fr"

try:
    fragment_length = int(config["fragment_length"])
except:
    fragment_length = 200

# IMPORTANT: When using snakemake with argument --config key=True, the
# string "True" is assigned to variable "key". Assigning a boolean value
# does not seem to be possible. Therefore, --config key=False will also
# return the boolean value True as bool("False") gives True.
# In contrast, within a configuration file config.yaml, assigment of boolean
# values is possible.
try:
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
except:
    trim = False
    fastq_dir = "FASTQ"

# IMPORTANT: When using snakemake with argument --config key=True, the
# string "True" is assigned to variable "key". Assigning a boolean value
# does not seem to be possible. Therefore, --config key=False will also
# return the boolean value True as bool("False") gives True.
# In contrast, within a configuration file config.yaml, assigment of boolean
# values is possible.
try:
    if config["fastqc"]:
        fastqc = True
    else:
        fastqc = False
except:
    fastqc = False

try:
    trim_galore_opts = config["trim_galore_opts"]
except:
    trim_galore_opts = "--stringency 2"

# IMPORTANT: When using snakemake with argument --config key=True, the
# string "True" is assigned to variable "key". Assigning a boolean value
# does not seem to be possible. Therefore, --config key=False will also
# return the boolean value True as bool("False") gives True.
# In contrast, within a configuration file config.yaml, assigment of boolean
# values is possible.
try:
    if config["dedup"]:
        dedup = True
    else:
        dedup = False
except:
    dedup = False

# IMPORTANT: When using snakemake with argument --config key=True, the
# string "True" is assigned to variable "key". Assigning a boolean value
# does not seem to be possible. Therefore, --config key=False will also
# return the boolean value True as bool("False") gives True.
# In contrast, within a configuration file config.yaml, assigment of boolean
# values is possible.
try:
    if config["properpairs"]:
        properpairs = True
    else:
        properpairs = False
except:
    properpairs = False

try:
    mapq = int(config["mapq"])
except:
    mapq = 0

try:
    bw_binsize = int(config["bw_binsize"])
except:
    bw_binsize = 10


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
