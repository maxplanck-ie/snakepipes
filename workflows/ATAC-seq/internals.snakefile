import glob
import os
import subprocess

## Require configuration file (samples.yaml)
samples = [ os.path.basename(f) for f in glob.glob('filtered_bam/*.filtered.bam') ]
samples = [ os.path.basename(x).replace('.filtered.bam','') for x in samples ]

if sample_info and not os.path.isfile(sample_info):
    print("ERROR: Cannot find sample info file! ("+sample_info+")\n")
    exit(1)

if sample_info and not check_sample_info_header(sample_info):
    print("ERROR: Please use 'name' and 'condition' as column headers in sample info file! ("+sample_info+")\n")
    exit(1)
