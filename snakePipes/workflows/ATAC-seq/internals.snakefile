import glob
import os
import subprocess

## Require configuration file (samples.yaml)
samples = [ os.path.basename(f) for f in glob.glob('filtered_bam/*.filtered.bam') ]
samples = [ os.path.basename(x).replace('.filtered.bam','') for x in samples ]

if sampleSheet and not os.path.isfile(sampleSheet):
    print("ERROR: Cannot find sample sheet file! ("+sampleSheet+")\n")
    exit(1)

if sampleSheet and not cf.check_sample_info_header(sampleSheet):
    print("ERROR: Please use 'name' and 'condition' as column headers in sample sheet file! ("+sampleSheet+")\n")
    exit(1)
