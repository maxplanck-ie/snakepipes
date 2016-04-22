import glob
import os

## basic variables #############################################################
try:
    indir = config["indir"]
except:
    indir = os.getcwd()

try:
    genome = config["genome"]
except:
    genome = None


## Setup table #################################################################

def get_samples_from_setup_table():
    """
    Parse setup table.
    """
    control = None
    samples_control = []
    samples_treatment = []
    with open(os.path.join(indir, "setup_table", "setup_table.tsv")) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            e = line.split()
            if not control:
                control = e[1]
            if e[1] == control:
                samples_control.append(e[0])
            else:
                samples_treatment.append(e[0])
    return(samples_control, samples_treatment)


samples_control, samples_treatment = get_samples_from_setup_table()
samples = samples_control + samples_treatment
# print(samples_control)
# print(samples_treatment)
# print(samples)


################################################################################
# print("Indir:", indir)
# print("Infiles:", infiles)
# print("Samples:", samples)
# print("Downsample:", downsample)
