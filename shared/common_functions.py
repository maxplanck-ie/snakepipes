### functions shared across workflows ##########################################
################################################################################

def get_fragment_length(infile):
    """
    Return median insert size from a metrics file created by
    Picard CollectInsertSizeMetrics
    Read 2 column text file, grep line by 1st column, return 2nd
    """
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("MEDIAN_INSERT_SIZE"):
                try:
                    median = next(f).split()[0]
                    return int(median)
                except:
                    print("ERROR: File", infile, "is NOT a proper Picard CollectInsertSizeMetrics metrics file.\n")
                    exit(1)
    # no match in infile
    print("ERROR: File", infile, "is NOT a proper Picard CollectInsertSizeMetrics metrics file.")
    exit(1)


########### Temp dir setup ###############

if "tempdir" in config:
    temp_path = config["tempdir"]+"/"

try: 
    output = subprocess.check_output("mktemp -d -p "+temp_path+"/ tmp.snakemake.XXXXXXXX",shell=True,stderr=subprocess.STDOUT)
    temp_path = output.decode().rstrip()+"/";
except subprocess.CalledProcessError:
    print("Failed to create temp dir under default temp path prefix ("+temp_path+")! Use "+outdir+" instead!")
    temp_path = outdir+"/"
