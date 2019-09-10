import glob
import os
import sys

## Main variables ##############################################################


### Functions ##################################################################

def readSampleSheet(sampleSheet, reads):
    """
    Read in preprocessing sample sheet and return a dictionary with keys
    the final sample prefixes and values the samples that should be merged
    """
    f = open(sampleSheet)
    pairedEnd = False
    ss = dict()
    for line in f:
        cols = line.strip().split("\t")
        if cols[2] not in ss:
            ss[cols[2]] = [[], []]
        if cols[1] == reads[0]:
            ss[cols[2]][0].append(cols[0])
        elif cols[1] == reads[1]:
            pairedEnd = True
            ss[cols[2]][1].append(cols[0])
        else:
            sys.exit("ERROR: The read mate designator for sample {} ({}) did not match either of those specified by --reads ({} and {})!\n".format(cols[2], cols[1], reads[0], reads[1]))

    # sanity check
    for k, v in ss.items():
        if len(v[1]) > 0:
            if len(v[0]) != len(v[1]):
                sys.exit("ERROR: There is a mismatch between the number of read 1 and read 2 entries in {}: [{}] and [{}]\n".format(k, v[0], v[1]))

    return ss, pairedEnd


### Initialization #############################################################
infiles = sorted(glob.glob(os.path.join(str(indir or ''), '*{}'.format(ext))))
samples = cf.get_sample_names(infiles, ext, reads)
pairedEnd = cf.is_paired(infiles, ext, reads)

if not samples and (not sampleSheet or sampleSheet == ""):
    sys.exit("Error! NO samples found in dir {}!!!\n".format(str(indir or '')))
