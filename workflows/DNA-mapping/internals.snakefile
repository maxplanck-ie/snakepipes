import glob
import os
import subprocess

## Main variables ##############################################################


### Functions ##################################################################

# When modifying the function update_filter(), double-check wether the rule
# samtools_filter has to be modified concordantly
def update_filter(samples,dedup,properpairs,mapq):
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

# ensure that only the specified filters are applied to all files
# delete already filtered BAM files if they were generated with different
# filtering parameters in previous runs
update_filter(samples,dedup,properpairs,mapq)
