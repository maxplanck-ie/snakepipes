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


### Variable defaults ##########################################################

try:
    indir = config["indir"]
except:
    indir = os.getcwd()

try:
    outdir = config["outdir"]
except:
    outdir = os.getcwd()

try:
    genome = config["genome"]
except:
    genome = None

try:
    bw_binsize = int(config["bw_binsize"])
except:
    bw_binsize = 10

try:
    fragment_length = int(config["fragment_length"])
except:
    fragment_length = 200

try:
    paired = config["paired"]
except:
    paired = False


### Initialization #############################################################

# TODO: catch exception if ChIP-seq samples are not unique
# read ChIP-seq dictionary from config.yaml:
# { ChIP1: { control: Input1, broad: True }, ChIP2: { control: Input2, broad: false }
chip_dict = config["chip_dict"]

# create unique sets of control samples, ChIP samples with and without control
control_samples = set()
chip_samples_w_ctrl = set()
chip_samples_wo_ctrl = set()
for chip_sample, value in chip_dict.items():
    if value['control']:
        control_samples.add(value['control'])
        chip_samples_w_ctrl.add(chip_sample)
    else:
        chip_samples_wo_ctrl.add(chip_sample)
    # set broad to False if not specified or set to False
    if not value['broad']:
        chip_dict[chip_sample]['broad'] = False

control_samples = list(sorted(control_samples))
chip_samples_w_ctrl = list(sorted(chip_samples_w_ctrl))
chip_samples_wo_ctrl = list(sorted(chip_samples_wo_ctrl))
chip_samples = sorted(chip_samples_w_ctrl + chip_samples_wo_ctrl)
all_samples = sorted(control_samples + chip_samples)

# TODO: do consistency check whether BAM files and .bam.bai index files of all samples given in chip_dict exist
