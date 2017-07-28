import glob
import os
import subprocess
import re
import yaml


## Main variables ##############################################################


### Functions ##################################################################

### Variable defaults ##########################################################

### Initialization #############################################################

# TODO: catch exception if ChIP-seq samples are not unique
# read ChIP-seq dictionary from config.yaml:
# { ChIP1: { control: Input1, broad: True }, ChIP2: { control: Input2, broad: false }
#config["atac_dict"] = {}

## Require configuration file (samples.yaml)
if not os.path.isfile(samples_config):
    print("ERROR: Cannot find samples file ("+samples_config+")")
    exit(1)

atac_dict = {}
with open(samples_config, "r") as f:
    atac_dict_tmp = yaml.load(f)
    if "atac_dict" in atac_dict_tmp and atac_dict_tmp["atac_dict"] :
        atac_dict = atac_dict_tmp["atac_dict"]
    else:
        print("\n  Error! Sample config has empty or no 'atac_dict' entry! ("+config["samples_config"]+") !!!\n\n")
        exit(1)
    del atac_dict_tmp

cf.write_configfile(os.path.join("atac_samples.yaml"),atac_dict)

# create unique sets of control samples, ChIP samples with and without control


samples = [x for x in [atac_dict]]



# consistency check whether all required files exist for all samples
for sample in samples:
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
