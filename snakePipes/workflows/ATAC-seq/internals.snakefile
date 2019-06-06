import glob
import os
import subprocess

## Require configuration file (samples.yaml)
samples = [ os.path.basename(f) for f in glob.glob('filtered_bam/*.filtered.bam') ]
samples = [ os.path.basename(x).replace('.filtered.bam','') for x in samples ]

if sampleSheet:
    cf.check_sample_info_header(sampleSheet)
    if not cf.check_replicates(sampleSheet):
        print("\nWarning! CSAW cannot be invoked without replicates!\n")
        sys.exit()

if not fromBAM:
    if paired:
        if not os.path.isfile(os.path.join(workingdir, "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv")):
            sys.exit('ERROR: {} is required but not present\n'.format(os.path.join(workingdir, "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv")))

    # consistency check whether all required files exist for all samples
    for sample in samples:
        req_files = [
            os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam"),
            os.path.join(workingdir, "filtered_bam/"+sample+".filtered.bam.bai")
            ]

        # check for all samples whether all required files exist
        for file in req_files:
            if not os.path.isfile(file):
                print('ERROR: Required file "{}" for sample "{}" specified in '
                      'configuration file is NOT available.'.format(file, sample))
                exit(1)

        
else:
    bamFiles = sorted(glob.glob(os.path.join(str(fromBAM or ''), '*' + bamExt)))
    bamSamples = cf.get_sample_names_bam(bamFiles, bamExt)
    
    
    bamDict = dict.fromkeys(bamSamples)
    
    print(bamFiles)
    print(bamSamples)
    
    mapping_prg="EXTERNAL_BAM"
    indir = fromBAM
    samples=bamSamples
    downsample = None
