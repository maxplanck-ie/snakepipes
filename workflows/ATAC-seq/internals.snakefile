import glob
import os
import subprocess
import re
import yaml


## Require configuration file (samples.yaml)
samples = [ os.path.basename(f) for f in glob.glob('filtered_bam/*.filtered.bam') ]
samples = [ os.path.basename(x).replace('.filtered.bam','') for x in samples ]
#print(samples)

# consistency check whether all required files exist for all samples
for sample in samples:
    req_files = [
        os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam"),
        os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam.bai")
    ]
    for file in req_files:
        if not os.path.isfile(file):
            print('ERROR: Required file "{}" for sample "{}" specified in '
                  'configuration file is NOT available.'.format(file, sample))
            exit(1)

if not os.path.exists(outdir_MACS2):
    os.makedirs(outdir_MACS2)

## do not remove X from norm?
ignoreForPeaks = ignore_forNorm.split(' ')
if "X" in ignoreForPeaks:
    ignoreForPeaks.remove("X")
ignoreForPeaksFile = os.path.join(outdir_MACS2,'Annotation_'+genome+"_chroms.grep_ignore")

with open(ignoreForPeaksFile,'w') as f:
    for chrom in ignoreForPeaks:
        f.write('^'+chrom +'\n')
