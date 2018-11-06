.. _HiC:

HiC
===

What it does
------------

The snakePipes HiC workflow allows users to process their HiC data from raw fastq files to
corrected HiC matrices and TADs. The workflow utilized mapping by BWA, followed by analysis
using `HiCExplorer` (https://www.nature.com/articles/s41467-017-02525-w). The output matrices
are currently in HiCExplorer .hdf5 format.

.. image:: ../images/HiC_pipeline.png

Input requirements and outputs
------------------------------

This pipeline requires paired-end reads fastq files as input in order to build a contact matrix and to call TADs.
Prior to building the matrix, the pipeline maps reads against a user-specified reference genome.
The output of mapping step is then used for building the contact matrix.

Workflow configuration file
---------------------------

Default parameters from the provided config file can be altered by user. Below is
the config file description for the HiC workflow :

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
      configfile:
      cluster_configfile:
      #set to true if running locally
      local: False
      #number of threads
      max_jobs: 5
      ## directory with fastq files
      indir:
      ## preconfigured target genomes (mm9,mm10,dm3,...) , see /path/to/snakemake_workflows/shared/organisms/
      ## Value can be also path to your own genome config file!
      genome:
      ## FASTQ file extension (default: ".fastq.gz")
      ext: '.fastq.gz'
      ## paired-end read name extension (default: ['_R1', "_R2"])
      reads: ['_R1', '_R2']
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
      # which restriction enzyme was used
      enzyme: HindIII
      # bin size in base pairs, if RF resolution is not required
      bin_size: 10000
      # build matrix only for a given region chr:start-end
      restrict_region:
      # Create files after merging a given number bins to be merged (default 0 = bins are not merged)
      nbins_toMerge: 0
      # shall we merge the samples?
      merge_samples: false
      # parameters for hicFindTADs
      tadparams: '--thresholdComparisons 0.01'
      # Should hicPlotDistVsCounts be run?
      distVsCount: false
      #Parameters to run hicPlotDistVsCounts
      distVsCountParams:
      #Terminate the pipeline before calling TADs
      noTAD: false
      #Terminate the pipeline before correting the matrices with a certain cutoff value
      noCorrect: false
      #Chromosomes of interest to build matrix on them
      chromosomes:
      # a .tsv file contains names and replicates of samples. It is needed if merge_samples
      sample_info:
      ################################################################################
      # Call snakemake directly, i.e. without using the wrapper script:
      #
      # Please save a copy of this config yaml file and provide an adjusted config
      # via '--configfile' parameter!
      # example call:
      #
      # snakemake --snakefile /path/to/snakemake_workflows/workflows/HiC/Snakefile
      #           --configfile /path/to/snakemake_workflows/workflows/HiC/defaults.yaml
      #           --directory /path/to/outputdir
      #           --cores 32
      ################################################################################


Structure of output directory
-----------------------------

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

- BWA folder contains the mapping results in `bam` format. The files were obtained after running BWA on each of the paired-end reads individually.
- FASTQ folder has soft links to the original input reads as well as a copy of trimmed reads if trimming has been called.
- HiC_matrices folder accommodates the contact matrices generated by hicBuildMatrix (https://www.nature.com/articles/s41467-017-02525-w). In case of merging samples or merging bins the initial matrix is saved in this folder along with the merged ones.
  - QCplot includes the QC measurements for each sample along with a diagnostic plot which illustrates a distribution of counts per bin.  This information can be used to set a cutoff to prune (correct) the contact matrix.
  The cutoff value is computed by the pipeline and by default will be applied to build a corrected matrix. Generated matrices by the pipeline can further be used for downstream analysis such as detecting A/B compartments and they can also be visualized using hicPlotMatrix.
- HiC_matrices_corrected folder is in fact containing the corrected matrix which has been generated after pruning as has been mentioned above.
- TADs folder includes the output of calling TADs using hicFindTADs (https://www.nature.com/articles/s41467-017-02525-w). The output contains TAD boundaries,  TAD domains and TAD scores. These along with the matrices can be visualized together as several tracks using pyGenomeTrack or can be dynamically viewed via hicBrowser.
  Check figure below as an example.
    
.. image:: ../images/HiC_tracks.png

Command line options
--------------------

.. argparse::
    :func: parse_args
    :filename: ../snakePipes/workflows/HiC/HiC
    :prog: HiC
    :nodefault:
