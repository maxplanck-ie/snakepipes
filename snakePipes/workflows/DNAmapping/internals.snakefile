import glob
import os
import subprocess
import warnings

## Main variables ##############################################################


### Functions ##################################################################


### Variable defaults ##########################################################
mode = list(map( str.strip, re.split(',|;', config["mode"]) ))
mode = [element.lower() for element in mode]
## genome names for allele-sp mapping
strains = list(map( str.strip, re.split(',|;', config["strains"]) ))
## trim
fastq_dir = "FASTQ"
if trim:
    fastq_indir_trim = "FASTQ"
    if trimmer == "trimgalore":
        fastq_dir = "FASTQ_TrimGalore"
    elif trimmer == "cutadapt":
        fastq_dir = "FASTQ_Cutadapt"
    elif trimmer == "fastp":
        fastq_dir = "FASTQ_fastp"


### Initialization #############################################################

if "allelic-whatshap" in mode:
    if "allelic-mapping" in mode:
        sys.exit("Allelic-mapping and allelic-whatshap modes are not compatible. Please choose one or another.")
    if not pvcf:
        sys.exit("Allelic-whatshap mode was specified but no phased vcf file was provided. Please provided the path to a phased vcf file.")
    if not os.path.isfile(pvcf):
        sys.exit(f"File {pvcf} doesn't exist.")



infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*'+ext)))
if infiles == []:
    sys.exit("Error! Samples extnesion in {} are not {}. "
             "Please change the extensions to it or update the config.yaml file "
             "with your desired extension.".format(indir,ext))
samples = cf.get_sample_names(infiles,ext,reads)

pairedEnd = cf.is_paired(infiles,ext,reads)

del infiles

if not samples:
    sys.exit("\n  Error! NO samples found in dir "+str(indir or '')+"!!!\n\n")

fromBAM = None

idxRange = 1
if pairedEnd:
    idxRange = 2

# clean up filtered_bam if needed appropriate
os.makedirs(outdir, exist_ok=True)
filt = ""
if dedup:
    filt += "-F 1024 "
    assert UMIDedup is False, "\nPlease use either --UMIDedup (UMI-tools dedup) or --dedup (via sambamba and samtools)!\n"
    "should be called!"
if properPairs:
    filt += "-f 2 "
if mapq is not None and mapq > 0:
    filt += "-q {} ".format(mapq)
filter_rules = os.path.join(outdir, "filter_rules")
if os.path.exists(filter_rules):
    f = open(filter_rules)
    cont = f.read()
    f.close()
    if cont != filt:
        f = open(filter_rules, "w")
        f.write(filt)
        f.close()
else:
    f = open(filter_rules, "w")
    f.write(filt)
    f.close()


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
        print("\n useSpikeInForNorm was specified but no spikein genome was detected - no spikeIn chromosomes found with extention " + spikeinExt + " .\n\n")
        exit(1)
