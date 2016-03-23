import glob
import os


## Functions ###################################################################

def get_sample_names(infiles):
    """
    Get sample names without extensions
    """
    s = []
    for x in sorted(glob.glob(os.path.join(indir, '*'+ext))):
        x = os.path.basename(x).replace(ext,"")
        try:
            x = x.replace(reads[0],"").replace(reads[1],"")
        except:
            pass
        s.append(x)
    return(sorted(list(set(s))))


def is_paired(infiles):
    """
    Check for paired input files
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
    return(paired)


def get_from_file(a, infile):
    """
    Read 2 column text file, grep line by 1st column, return 2nd
    """
    b = None
    with open(infile,"r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(a):
                try:
                    b = line.split()[1]
                except:
                    b = ""
                break
    return(b)


## Variable defaults ###########################################################

try:
    indir = config["indir"]
except:
    indir = os.getcwd()

## paired reads extensions
try:
    reads = config["reads"]
except:
    reads = ["_R1", "_R2"]

## FASTQ file extension
try:
    ext = config["ext"]
except:
    ext = ".fastq.gz"

## Downsample (number of reads)
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

try:
    fragment_length = int(config["fragment_length"])
except:
    fragment_length = 200

try:
    if config["trim"] == "True":
        trim = True
except:
    trim = False

try:
    trim_galore_opts = config["trim_galore_opts"]
except:
    trim_galore_opts = "--stringency 2"

try:
    binsize = int(config["binsize"])
except:
    binsize = 25

## Default input directory for the mapping program
if trim:
    fastq_dir = "TrimGalore"
else:
    fastq_dir = "FASTQ"


infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = get_sample_names(infiles)

paired = is_paired(infiles)

if not paired:
    reads = [""]


################################################################################
# print("Indir:", indir)
# print("Infiles:", infiles)
# print("Samples:", samples)
# print("Downsample:", downsample)
