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

if not fromBam:
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
    bamFiles = sorted(glob.glob(os.path.join(str(fromBam or ''), '*'+bam_ext)))
    bamSamples = cf.get_sample_names_bam(bamFiles,bam_ext)
    
    
    bamDict = dict.fromkeys(bamSamples)
    
    print(bamFiles)
    print(bamSamples)
    
    mapping_prg="EXTERNAL_BAM"
    indir = fromBam
    samples=bamSamples
    downsample = None