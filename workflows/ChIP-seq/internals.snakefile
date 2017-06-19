import glob
import os
import subprocess
import re
import yaml


## Main variables ##############################################################


### Functions ##################################################################

def get_control(sample):
    """
    Return control sample name for a given ChIP-seq sample
    Return False if given ChIP-seq sample has no control
    """
    if sample in chip_samples_w_ctrl:
        return chip_dict[sample]['control']
    else:
        return False


def is_broad(sample):
    """
    Return True if given ChIP-seq sample is annotated as sample with
    broad enrichment, else return False
    """
    if sample in chip_dict:
        return chip_dict[sample]['broad']
    else:
        return False


def is_chip(sample):
    """
    Return True if a given sample is a ChIP-seq sample
    Else return False
    """
    return (sample in chip_samples)


### Variable defaults ##########################################################
### Initialization #############################################################

# TODO: catch exception if ChIP-seq samples are not unique
# read ChIP-seq dictionary from config.yaml:
# { ChIP1: { control: Input1, broad: True }, ChIP2: { control: Input2, broad: false }

with open(config["samples_config"], "r") as f:
    config["chip_dict"] = yaml.load(f)["chip_dict"]

chip_dict = config["chip_dict"]
cf.write_configfile(os.path.join(outdir,"chip_samples.yaml"),chip_dict)

# create unique sets of control samples, ChIP samples with and without control
control_samples = set()
chip_samples_w_ctrl = set()
chip_samples_wo_ctrl = set()
for chip_sample, value in chip_dict.items():
    # set control to False if not specified or set to False
    if 'control' not in chip_dict[chip_sample] or not value['control']:
        chip_dict[chip_sample]['control'] = False
        chip_samples_wo_ctrl.add(chip_sample)
    else:
        control_samples.add(value['control'])
        chip_samples_w_ctrl.add(chip_sample)
    # set broad to False if not specified or set to False
    if 'broad' not in chip_dict[chip_sample] or not value['broad']:
        chip_dict[chip_sample]['broad'] = False

control_samples = list(sorted(control_samples))
chip_samples_w_ctrl = list(sorted(chip_samples_w_ctrl))
chip_samples_wo_ctrl = list(sorted(chip_samples_wo_ctrl))
chip_samples = sorted(chip_samples_w_ctrl + chip_samples_wo_ctrl)
all_samples = sorted(control_samples + chip_samples)

# consistency check whether all required files exist for all samples
for sample in all_samples:
    req_files = [
        os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam"),
        os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam.bai"),
        os.path.join(workingdir, "Picard_qc/MarkDuplicates/"+sample+".mark_duplicates_metrics.txt"),
        os.path.join(workingdir, "Picard_qc/AlignmentSummaryMetrics/"+sample+".alignment_summary_metrics.txt")
    ]
    if paired:
        req_files.append(os.path.join(workingdir, "Picard_qc/InsertSizeMetrics/"+sample+".insert_size_metrics.txt"))

    # check for all samples whether all required files exist
    for file in req_files:
        if not os.path.isfile(file):
            print('ERROR: Required file "{}" for sample "{}" specified in '
                  'configuration file is NOT available.'.format(file, sample))
            exit(1)
