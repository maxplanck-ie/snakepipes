### functions shared across workflows ##########################################
################################################################################
import subprocess

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


def make_temp_dir(tempdir, fallback_dir, verbose=False):
    try:
        output = subprocess.check_output("mktemp -d -p "+tempdir+"/ tmp.snakemake.XXXXXXXX",shell=True,stderr=subprocess.STDOUT)
        temp_path = output.decode().rstrip()+"/";
    except subprocess.CalledProcessError:
        try:
            print("\nFailed to create temp dir under temp path prefix ("+tempdir+")! Try fallback: "+fallback_dir+"/ ...")
            output = subprocess.check_output("mktemp -d -p "+fallback_dir+"/ tmp.snakemake.XXXXXXXX",shell=True,stderr=subprocess.STDOUT)
            temp_path = output.decode().rstrip()+"/";
        except subprocess.CalledProcessError:
            print("\nAlso failed to create temp dir under fallback prefix ("+fallback_dir+"/)!")
            exit(1)
    if verbose:
        print("\ntemp dir created: "+temp_path)
    return temp_path
