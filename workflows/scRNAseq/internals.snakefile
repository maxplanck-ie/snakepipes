import glob
import os
import subprocess


## Main variables ##############################################################


### Functions ###################################################################


### Variable defaults ##########################################################

## set trimming related dirs correctly
## we trim only R2 with transfered barcode info under "FASTQ_barcoded"
## only cutadapt is supported right now due to poly-A trimming 

fastq_dir = "FASTQ_barcoded"

if trim:
	fastq_indir_trim = "FASTQ_barcoded"
	if trim_prg == "trimgalore":
		fastq_dir = "FASTQ_TrimGalore"
	elif trim_prg == "cutadapt":
		fastq_dir = "FASTQ_Cutadapt"

mapping_prg = "STAR_genomic"

star_options = "--twopassMode Basic"
hisat_options = ""

### Initialization #############################################################

infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = cf.get_sample_names(infiles,ext,reads)

## we just check if we have correctly paired fastq files
if not cf.is_paired(infiles,ext,reads):
    print("This workflow requires paired-end read data!")
    exit(1)

## After barcode transfer to R2 we have only single end data / R2
## but we need to keep "reads" for rule fastq_barcode
paired = False

### barcode pattern extraction #################################################
pattern = re.compile("[N]+")

if pattern.search(barcode_pattern) is not None:
	UMI_offset = pattern.search(barcode_pattern).start() + 1 
	UMI_length = pattern.search(barcode_pattern).end() - UMI_offset + 1
else:
	print("Provided barcode pattern does not contain any 'N'! Exit...\n")
	exit(1)

pattern = re.compile("[X]+")

if pattern.search(barcode_pattern) is not None:
	CELLI_offset = pattern.search(barcode_pattern).start() + 1
	CELLI_length = pattern.search(barcode_pattern).end() - CELLI_offset + 1 
else:
	print("Provided barcode pattern does not contain any 'X'! Exit...\n")
	exit(1)
