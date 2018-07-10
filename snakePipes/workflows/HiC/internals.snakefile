import glob
import os
import subprocess
import pandas as pd
## Main variables ##############################################################


# minimum distance to consider between restriction sites. RS that are closer than this
# distance are merged. This value should be similar to the read length
MIN_RS_DISTANCE = 150
# Maximum distance of that a read can be from a restriction site. Reads farther than this distance are discarded.
# This number should be related to the higher end of fragment length distribution.
MAX_RS_DISTANCE = 1000
# this is to identify so called dangling ends. Option are 'AAGCTT' or 'GATC'
#fastq_dir = 'fastq'
#samples = ['SRR2240738', 'SRR2240739', 'SRR2240740']

### Functions ##################################################################

## find out which resolution is it (RF or binsize)
def get_matrixFile_suffix():
    if(RF_resolution is True):
        return("rf")
    else:
        return("bs")
matrixFile_suffix = get_matrixFile_suffix()

## seq dict for restriction enzymes
seq_dict = {
           'DpnII' : ['GATC', 'GATC'],
           'HindIII' : ['AAGCTT', 'AGCTT']
           }

## get sequences for restriction enzymes
def get_restriction_seq(name):
    seq = seq_dict[name][0]
    return(seq)

def get_dangling_seq(name):
    seq = seq_dict[name][1]
    return(seq)

## merge all the HiC matrices and create new file if asked
def merge_inputs():
    if(merge_samples is True):
        return(["HiC_matrices/all_matrices_merged.h5"#,"HiC_matrices/QCplots/all_matrices_merged_diagnostic_plot.pdf"
               ])
    else:
        return([])

## merge hic Bins if asked and create new files
def get_merged_bins():
    if(nbins_toMerge != 0):
        return(expand("HiC_matrices/{sample}_"+matrixFile_suffix+"_m"+str(nbins_toMerge)+".h5", sample=samples))
    else:
        return([])

## get MAD thresholds for the matrix correction
def get_mad_score(madfile):
    with open(madfile) as md:
        lower = 0.0
        for line in md:
            lower = float(line.split()[2])
    upper = -(3*lower)
    cutoff = str(lower) + " " + str(upper)
    return(cutoff)

## get sample grouping information
def get_sampleInfo(sample_info):
    sample_dict = {}
    if sample_info:     #Read the sample info and make a dictionary
        sample_conditions =  pd.read_csv(os.path.join(os.path.abspath(sample_info)), sep = '\t', header = 0)
        for id, row in sample_conditions.iterrows():
            v, k = row
            if k in sample_dict:
                sample_dict[k] = [sample_dict[k], v]
            else:
                sample_dict[k]=v
    else:
        for sample in samples:
            sample_dict['cond1'] = sample

    return sample_dict

## trim
fastq_dir = "FASTQ"
if trim:
    fastq_indir_trim = "FASTQ"
    if trim_prg == "trimgalore":
        fastq_dir = "FASTQ_TrimGalore"
    elif trim_prg == "cutadapt":
        fastq_dir = "FASTQ_Cutadapt"


### Initialization #############################################################
infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*'+ext)))
samples = cf.get_sample_names(infiles,ext,reads)
paired = cf.is_paired(infiles,ext,reads)
del infiles

if not samples:
    print("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")
    exit(1)

if not paired:
    print("\n  Error! Paired-end samples not detected. "
          "Hi-C workflow requires paired-end samples "+str(indir or '')+"!!!\n\n")
    exit(1)
