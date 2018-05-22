import glob
import os


## Require configuration file (samples.yaml)
samples = [ os.path.basename(f) for f in glob.glob('filtered_bam/*.filtered.bam') ]
samples = [ os.path.basename(x).replace('.filtered.bam','') for x in samples ]
