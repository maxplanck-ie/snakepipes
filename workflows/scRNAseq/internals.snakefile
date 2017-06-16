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

#try:
#	filter_annotation = config["filter_annotation"]
#except:
#    filter_annotation = "''"
	

### Initialization #############################################################

infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = cf.get_sample_names(infiles,ext,reads)

## we just check if we have correctly paired fastq files
if not cf.is_paired(infiles,ext,reads):
    print("This workflow requires paired-end read data!")
    exit(1)

## After barcode transfer to R2 we have only single end data / R2
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

print("UMI_LEN:",UMI_length,"  UMI_offset:",UMI_offset,"\n")
print("CELLI_LEN:",CELLI_length,"  CELLI_offset:",CELLI_offset,"\n")
