import glob
import os
import subprocess
import re
import yaml
import sys
import pandas as pd
import warnings

### Functions ##################################################################

def get_control(sample):
    """
    Return control sample name for a given ChIPseq sample
    Return False if given ChIPseq sample has no control
    """
    if sample in chip_samples_w_ctrl:
        return chip_dict[sample]['control']
    else:
        return False


def get_control_name(sample):
    """
    Return control sample alias for a given ChIPseq sample
    Return False if given ChIPseq sample has no control
    """
    if sample in chip_samples_w_ctrl:
        if 'control' in chip_dict[sample] and chip_dict[sample]['control'] != None:
            return chip_dict[sample]['control']
        else:
            return False
    else:
        return False


def is_broad(sample):
    """
    Return True if given ChIPseq sample is annotated as sample with
    broad enrichment, else return False
    """
    if sample in chip_dict:
        return chip_dict[sample]['broad']
    else:
        return False


def is_chip(sample):
    """
    Return True if a given sample is a ChIPseq sample
    Else return False
    """
    return (sample in chip_samples)


def is_allelic(workingdir):
    if os.path.isdir(os.path.join(workingdir,'allelic_bams') ) and os.listdir(os.path.join(workingdir,'allelic_bams') ) != []:
        return True
    else:
        return False


def get_pe_frag_length(sample, frag_len_file):
    try:
        df = pd.read_csv(frag_len_file, header = None, skiprows = 1, sep = "\t")
        df = df.loc[df[0] == sample]
        frag_len = int(df[5].values[0])
        return str(frag_len)
    except:
        warnings.warn("fragmentSize.metric.tsv is empty, this sets "
                      "--extsize of MACS2 to an empty string. Fix this and run MACS2 again!")
        return " "
### Variable defaults ##########################################################
### Initialization #############################################################

allele_info=is_allelic(workingdir)

# TODO: catch exception if ChIPseq samples are not unique
# read ChIPseq dictionary from config.yaml:
# { ChIP1: { control: Input1, broad: True }, ChIP2: { control: Input2, broad: false }
#config["chip_dict"] = {}

if not os.path.isfile(samples_config):
    print("ERROR: Cannot find samples file ("+samples_config+")")
    exit(1)

if sampleSheet:
    cf.check_sample_info_header(sampleSheet)
    if not cf.check_replicates(sampleSheet):
        print("\nWarning! CSAW cannot be invoked without replicates and will not be run!\n")
        if not peakCaller=="Genrich":
            sys.exit()
    isMultipleComparison = cf.isMultipleComparison(sampleSheet)
else:
    isMultipleComparison = False

chip_dict = {}
with open(samples_config, "r") as f:
    chip_dict_tmp = yaml.load(f, Loader=yaml.FullLoader)
    if "chip_dict" in chip_dict_tmp and chip_dict_tmp["chip_dict"] :
        chip_dict = chip_dict_tmp["chip_dict"]
    else:
        print("\n  Error! Sample config has empty or no 'chip_dict' entry! ("+config["samples_config"]+") !!!\n\n")
        exit(1)
    del chip_dict_tmp

cf.write_configfile(os.path.join("chip_samples.yaml"), chip_dict)

# create unique sets of control samples, ChIP samples with and without control
control_samples = set()
chip_samples_w_ctrl = set()
chip_samples_wo_ctrl = set()
for chip_sample, value in chip_dict.items():
    # set control to False if not specified or set to False
    if 'control' not in chip_dict[chip_sample] or value['control'] is None:
        chip_dict[chip_sample]['control'] = False
        chip_samples_wo_ctrl.add(chip_sample)
    else:
        control_samples.add(value['control'])
        chip_samples_w_ctrl.add(chip_sample)
    # set broad to False if not specified or set to False
    if 'broad' not in chip_dict[chip_sample] or not value['broad']:
        chip_dict[chip_sample]['broad'] = False


control_samples = list(sorted(control_samples))
# get a list of corresp control_names for chip samples
control_names = []
for chip_sample in chip_samples_w_ctrl:
    control_names.append(get_control_name(chip_sample))

chip_samples_w_ctrl = list(sorted(chip_samples_w_ctrl))
chip_samples_w_ctrl = list(filter(None, chip_samples_w_ctrl))
chip_samples_wo_ctrl = list(sorted(chip_samples_wo_ctrl))
chip_samples_wo_ctrl = list(filter(None, chip_samples_wo_ctrl))
chip_samples = sorted(chip_samples_w_ctrl + chip_samples_wo_ctrl)
chip_samples = list(filter(None, chip_samples))
all_samples = sorted(control_samples + chip_samples) if control_samples else chip_samples
all_samples = list(filter(None, all_samples))

#useful for debugging purposes; change the mode to info or debug
if chip_samples_wo_ctrl:
    warnings.warn( str(len(chip_samples_wo_ctrl)) + " out of " + str(len(chip_samples)) + " have no matching control ")
if chip_samples_w_ctrl:
    warnings.warn( str(len(chip_samples_w_ctrl)) + " out of " + str(len(chip_samples)) + " have a matching control ")

if not fromBAM:
    if pairedEnd and not useSpikeInForNorm:
        if not os.path.isfile(os.path.join(workingdir, "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv")):
            sys.exit('ERROR: {} is required but not present\n'.format(os.path.join(workingdir, "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv")))

    # consistency check whether all required files exist for all samples
    for sample in all_samples:
        req_files = [
            os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam"),
            os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam.bai"),
            ]
        if allele_info:
            req_files.append(os.path.join(workingdir, "bamCoverage/allele_specific/"+sample+".genome1.seq_depth_norm.bw"))
        else:
            if not useSpikeInForNorm:
                req_files.append(os.path.join(workingdir, "bamCoverage/"+sample+".filtered.seq_depth_norm.bw"))

        # check for all samples whether all required files exist
        for file in req_files:
            if not os.path.isfile(file):
                print('ERROR: Required file "{}" for sample "{}" specified in '
                      'configuration file is NOT available.'.format(file, sample))
                exit(1)


else:
    bamFiles = sorted(glob.glob(os.path.join(str(fromBAM or ''), '*' + bamExt)))
    bamSamples = cf.get_sample_names_bam(bamFiles, bamExt)

    bamDict = dict.fromkeys(bamSamples)

    for sample in all_samples:
        if sample not in bamDict:
            sys.exit("No bam file found for chip sample {}!".format(sample))
    aligner = "EXTERNAL_BAM"
    indir = fromBAM
    downsample = None

samples = all_samples
if not samples:
    print("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")
    exit(1)


##filter sample dictionary by the subset of samples listed in the 'name' column of the sample sheet
def filter_dict(sampleSheet,input_dict):
    f=open(sampleSheet,"r")
    nameCol = None
    nCols = None
    names_sub=[]
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            nameCol = cols.index("name")
            nCols = len(cols)
            continue
        elif idx == 1:
            if len(cols) - 1 == nCols:
                nameCol += 1
        if not len(line.strip()) == 0:
            names_sub.append(line.split('\t')[nameCol])
    f.close()
    output_dict = dict((k,v) for k,v in input_dict.items() if k in names_sub)
    return(output_dict)

if sampleSheet:
    if chip_samples_w_ctrl:
        filtered_dict = filter_dict(sampleSheet,dict(zip(chip_samples_w_ctrl, [ get_control_name(x) for x in chip_samples_w_ctrl ])))
    else:
        filtered_dict = filter_dict(sampleSheet,dict(zip(chip_samples_wo_ctrl, [None]*len(chip_samples_wo_ctrl))))
    genrichDict = cf.sampleSheetGroups(sampleSheet,isMultipleComparison)
    if not isMultipleComparison:
        for k in genrichDict.keys():
            genrichDict[k]=[item for item in genrichDict[k] if item in chip_samples]
        reordered_dict = {k: filtered_dict[k] for k in [item for sublist in genrichDict.values() for item in sublist]}
    else:
        #print(genrichDict)
        reordered_dict = {}
        for g in genrichDict.keys():
            for k in genrichDict[g].keys():
                genrichDict[g][k]=[item for item in genrichDict[g][k] if item in chip_samples]
                reordered_dict[g] = {k: filtered_dict[k] for k in [item for sublist in genrichDict[g].values() for item in sublist]}
else:
    genrichDict = {"all_samples": chip_samples}



#################### functions and checks for using a spiked-in genome for normalization ########################################
def check_if_spikein_genome(genome_index,spikeinExt):
    resl=[]
    if os.path.isfile(genome_index):
        with open(genome_index) as ifile:
            for line in ifile:
                resl.append(re.search(spikeinExt, line))
        if any(resl):
            warnings.warn("\n Spikein genome detected - at least one spikeIn chromosome found with extention " + spikeinExt + " .\n\n")
            return True
        else:
            return False
    else:
        print("\n  Error! Genome index file "+ genome_index +" not found!!!\n\n")
        exit(1)

def get_host_and_spikein_chromosomes(genome_index, spikeinEx):
    hostl=dict()
    spikeinl=dict()
    with open(genome_index) as ifile:
        for line in ifile:
            try:
                entry = line.split('\t')[0]
                length = line.split('\t')[1]
                if re.search(spikeinExt, entry):
                    spikeinl[entry] = length
                else:
                    hostl[entry] = length
            except:
                warnings.warn("check for empty lines in the index file!")
                continue
    return([hostl,spikeinl])

if useSpikeInForNorm:
    part=['host','spikein']
    spikein_detected=check_if_spikein_genome(genome_index,spikeinExt)
    if spikein_detected:
        host_chr, spikein_chr =get_host_and_spikein_chromosomes(genome_index,spikeinExt)
        spikein_region = ""
        if len(spikein_chr.items()) == 1:
            k, v = next(iter(spikein_chr.items()))
            spikein_region = ":0:".join([str(k),str(v)])
    else:
        print("\n No spikein genome detected - no spikeIn chromosomes found with extention " + spikeinExt + " .\n\n")
        exit(1)

if externalBed:
    if os.path.isfile(externalBed):
        peakCaller = os.path.splitext(os.path.basename(externalBed))[0]
    else:
        warnings.warn("{} file not found.".format(externalBed))
