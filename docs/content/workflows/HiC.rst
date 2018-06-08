.. _HiC:

HiC
============

Input requirements and outputs:
-------------------------------------------
This pipeline requires paired-end reads fastq files as input in order to build a contact matrix and to call TADs. 
Prior to building a matrix, the pipeline maps reads against a user provided reference genome. The output of mapping step is then used for building the contact matrix.

.. argparse::
   :filename: ../workflows/HiC/HiC
   :func: parse_args
   :prog: HiC

Usage example
------------------------------
.. parsed-literal::

   $ /path/to/snakemake/workflows/HiC/HiC [options] -i <input.dir> -o <output.dir>

Workflow config file
----------------------------------
Default parameters from the provided config file can be altered by user. Note that proper reads name extension from the fastq files should be entered as values of the 'reads' key of the config file. 

.. parsed-literal::

 ################################################################################
 # This file is the default configuration of the HiC workflow!
 #
 # In order to adjust some parameters, please either use the wrapper script
 # (eg. /path/to/snakemake_workflows/workflows/HiC/HiC)
 # or save a copy of this file, modify necessary parameters and then provide
 # this file to the wrapper or snakmake via '--configfile' option
 # (see below how to call the snakefile directly)
 #
 # Own parameters will be loaded during snakefile executiuon as well and hence
 # can be used in new/extended snakemake rules!
 ################################################################################
 ## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
 pipeline: hic
 outdir:
 onfigfile:
 luster_configfile: 
 local: False
 max_jobs: 5
 snakemake_options: '--use-conda --conda-prefix /path/to/conda/envs '
 tempdir: /data/extended/
 ## directory with fastq files
 indir:
 ## preconfigured target genomes (mm9,mm10,dm3,...) , see /path/to/snakemake/shared/organisms/
 ## Value can be also path to your own genome config file!
 genome:
 ## FASTQ file extension (default: ".fastq.gz")
 ext: .fastq.gz
 ## paired-end read name extension (default: ['_R1', "_R2"]) 
 reads: [_R1, _R2]
 mapping_prg: BWA
 ## Number of reads to downsample from each FASTQ file
 downsample:
 ## Options for trimming and fastqc
 trim: False
 trim_prg: cutadapt
 trim_options:
 fastqc: false
 verbose: False
 ## is the Matrix RF resolution?
 RF_resolution: false
 # which enzyme was used
 enzyme: HindIII
 # bin size in base pairs, if RF resolution is not required
 bin_size: 10000
 # build matrix only for given region
 restrict_region:
 # Create files with N merged bins (0 = No)
 nbins_toMerge: 0
 # shall we merge the samples?
 merge_samples: false
 # parameters for hicFindTADs
 tadparams: --thresholdComparisons 0.01
 distVsCount: false
 distVsCountParams: 
 noTAD: false
 noCorrect: false

 ################################################################################
 # Call snakemake directly, i.e. without using the wrapper script:
 #
 # Please save a copy of this config yaml file and provide an adjusted config
 # via '--configfile' parameter!
 # example call:
 #
 # snakemake --snakefile /path/to/snakemake/workflows/HiC/Snakefile
 #           --configfile /path/to/snakemake/workflows/HiC/defaults.yaml
 #           --directory /path/to/outputdir
 #           --cores 32
 ################################################################################

Structure of output directory
-------------------------------

The HiC pipeline will generate output of the following structure:

.. code:: bash

    $ tree -d -L 2 outputdir
  
output directory::

    outputdir
    |--BWA
    |--FASTQ
    |--HiC_matrices
    |   |--logs
    |   |--QCplots
    |--HiC_matrices_corrected
    |   |--logs
    |--TADs
        |--logs

