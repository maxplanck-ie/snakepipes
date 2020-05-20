import glob
import os
import subprocess
import sys


## Main variables ##############################################################


### Functions ###################################################################


### Variable defaults ##########################################################

## set trimming related dirs correctly
## we trim only R2 with transfered barcode info under "FASTQ_barcoded"
## only cutadapt is supported right now due to poly-A trimming 

if mode == "Gruen":
    fastq_dir = "FASTQ_barcoded"
    if trim:
        fastq_indir_trim = "FASTQ_barcoded"
        if trimmer == "trimgalore":
            fastq_dir = "FASTQ_TrimGalore"
        elif trimmer == "cutadapt":
            fastq_dir = "FASTQ_Cutadapt"
        elif trimmer == "fastp":
            fastq_dir = "FASTQ_fastp"
    else:
        fastq_indir_trim = None	
    aligner = "STAR_genomic"
elif mode == "STARsolo":
    fastq_indir_trim = None
    fastq_dir = "originalFASTQ"
    aligner = "STARsolo"

### Initialization #############################################################

infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = cf.get_sample_names(infiles,ext,reads)
fromBAM = None

if not samples:
    print("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")
    exit(1)

## we just check if we have correctly paired fastq files
if not cf.is_paired(infiles,ext,reads):
    print("This workflow requires paired-end read data!")
    exit(1)

## print deprecation message for modeGruen
if mode=="Gruen":
    print("Warning: mode Gruen is going to be deprecated by the end of 2020!")


if mode == "STARsolo":
    if myKit == "10Xv2":
        BCwhiteList = os.path.join(maindir,"workflows","scRNAseq","10x_737K-august-2016.txt")
        STARsoloCoords = ["17","10","1","16"]
    elif myKit == "10Xv3":
        BCwhiteList = "<( gzip -dc " + os.path.join(maindir,"workflows","scRNAseq","10x_3M-february-2018.txt.gz") + " )"
        STARsoloCoords = ["17","12","1","16"]
    elif myKit == "CellSeq192":
        BCwhiteList = os.path.join(maindir,"workflows","scRNAseq","celseq_barcodes.192.1col.txt")
        STARsoloCoords = ["1","6","7","6"]
    elif myKit == "CellSeq384":
        BCwhiteList = os.path.join(maindir,"workflows","scRNAseq","celseq_barcodes.384.1col.txt")
        STARsoloCoords = ["1","7","8","7"]
    elif myKit == "Custom":
        if not os.path.isfile(BCwhiteList):
            print("Provided barcode whitelist file doesn't exist! Exit...\n")
            exit(1)
  
## After barcode transfer to R2 we have only single end data / R2
## but we need to keep "reads" for rule fastq_barcode
pairedEnd = False
## we swap read extensions as we continue in SE mode but with R2
##some rules use a hardcoded reads[0] for SE
reads = reads[::-1]

### barcode pattern extraction #################################################
pattern = re.compile("[N]+")

if pattern.search(cellBarcodePattern) is not None:
    UMI_offset = pattern.search(cellBarcodePattern).start() + 1 
    UMI_length = pattern.search(cellBarcodePattern).end() - UMI_offset + 1
else:
    print("Provided barcode pattern does not contain any 'N'! Exit...\n")
    exit(1)

pattern = re.compile("[X]+")

if pattern.search(cellBarcodePattern) is not None:
    CELLI_offset = pattern.search(cellBarcodePattern).start() + 1
    CELLI_length = pattern.search(cellBarcodePattern).end() - CELLI_offset + 1 
else:
    print("Provided barcode pattern does not contain any 'X'! Exit...\n")
    exit(1)
