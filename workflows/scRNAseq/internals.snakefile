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

# IMPORTANT: When using snakemake with argument --config key=True, the
# string "True" is assigned to variable "key". Assigning a boolean value
# does not seem to be possible. Therefore, --config key=False will also
# return the boolean value True as bool("False") gives True.
# In contrast, within a configuration file config.yaml, assigment of boolean
# values is possible.
fastq_dir = ""
fastq_indir_trim = ""
trim = False

try:
    if config["trim"] == "trimgalore":
        trim = True
        fastq_dir = "FASTQ_TrimGalore"
        fastq_indir_trim = "FASTQ_barcoded"
    elif config["trim"] == "cutadapt":
        trim = True
        fastq_dir = "FASTQ_Cutadapt"
        fastq_indir_trim = "FASTQ_barcoded"
    else:
        trim = False
        fastq_dir = "FASTQ_barcoded"
except:
    trim = False
    fastq_dir = "FASTQ_barcoded"

try:
    trim_galore_opts = config["trim_galore_opts"]
except:
    trim_galore_opts = "-a \"A{{30}}\""


try:
    bw_binsize = int(config["bw_binsize"])
except:
    bw_binsize = 10

try:
    barcode_pattern = config["barcode_pattern"]
except:
    barcode_pattern = "NNNNNNXXXXXX"

try:
	barcode_file = config["barcode_file"]
except:
	barcode_file = workflow.basedir+"/celseq_barcodes.192.txt"
    
#try:
#	transcripts_exclude = config["transcripts_exclude"]
#except:
#	transcripts_exclude = "full|decay|miRNA|misc_RNA|snRNA|snoRNA|scaRNA|sRNA|scRNA|rRNA|pseudogene|3prime_overlapping_ncRNA|processed_transcript"

try:
	filter_annotation = config["filter_annotation"]
except:
    filter_annotation = "''"

    
### Initialization #############################################################

infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = get_sample_names(infiles)

if not is_paired(infiles):
    print("scRNA-seq needs currently paired read data!")
    exit(1)

## After barcode transfer to R2 we have only single end data / R2
paired = False

#paired = is_paired(infiles)

#if not paired:
#    reads = [""]

### barcode pattern extraction
pattern = re.compile("[N]+")
UMI_offset = pattern.search(barcode_pattern).start() + 1 
UMI_length = pattern.search(barcode_pattern).end() - UMI_offset + 1

pattern = re.compile("[X]+")
CELLI_offset = pattern.search(barcode_pattern).start() + 1
CELLI_length = pattern.search(barcode_pattern).end() - CELLI_offset + 1 


print("UMI_LEN:",UMI_length,"  UMI_offset:",UMI_offset,"\n")
print("CELLI_LEN:",CELLI_length,"  CELLI_offset:",CELLI_offset,"\n")
