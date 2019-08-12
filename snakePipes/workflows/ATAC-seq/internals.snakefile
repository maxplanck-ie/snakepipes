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

if fromBAM:
    bamFiles = sorted(glob.glob(os.path.join(str(fromBAM or ''), '*' + bamExt)))
    bamSamples = cf.get_sample_names_bam(bamFiles, bamExt)
    
    
    bamDict = dict.fromkeys(bamSamples)
    
    print(bamFiles)
    print(bamSamples)
    
    aligner = "EXTERNAL_BAM"
    indir = fromBAM
    samples = bamSamples
    downsample = None

    if pairedEnd:
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

##filter sample dictionary by the subset of samples listed in the 'name' column of the sample sheet
def filter_dict(sampleSheet):
    f=open(sampleSheet,"r")
    nameCol = None
    nCols = None
    names_sub=[]
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            nameCol = cols.index("name")
            nCols = len(cols)
            continue
        elif idx == 1:
            if len(cols) - 1 == nCols:
                nameCol += 1
        if not len(line.strip()) == 0:
            names_sub.append(line.split('\t')[nameCol])      
    f.close()
    output_dict = dict(zip(names_sub, [""]*len(names_sub)))
    return(output_dict)

if sampleSheet:
    filtered_dict = filter_dict(sampleSheet)
