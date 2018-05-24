.. _DNA-mapping:

DNA-mapping
============

Configuration file
------------------------

.. code:: bash

    $ cat snakemake_workflows/workflows/DNA-mapping/defaults.yaml

.. parsed-literal::

	################################################################################
	# This file is the default configuration of the DNA-mapping workflow!
	#
	# In order to adjust some parameters, please either use the wrapper script
	# (eg. /path/to/snakemake_workflows/workflows/DNA-mapping/DNA-mapping)
	# or save a copy of this file, modify necessary parameters and then provide
	# this file to the wrapper or snakmake via '--configfile' option
	# (see below how to call the snakefile directly)
	#
	# Own parameters will be loaded during snakefile executiuon as well and hence
	# can be used in new/extended snakemake rules!
	################################################################################
	## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
	outdir:
	configfile:
	local: False
	max_jobs: 5
	snakemake_options:
	tempdir: /data/extended/
	## directory with fastq files
	indir:
	## preconfigured target genomes (mm9,mm10,dm3,...) , see /path/to/snakemake_workflows/shared/organisms/
	## Value can be also path to your own genome config file!
	genome:
	## FASTQ file extension (default: ".fastq.gz")
	ext: .fastq.gz
	## paired-end read name extension (default: ['_R1', "_R2"])
	reads: [_R1, _R2]
	## Number of reads to downsample from each FASTQ file
	downsample:
	## Options for trimming
	trim: False
	trim_prg: cutadapt
	trim_options:
	## Bin size of output files in bigWig format
	bw_binsize: 25
	## Run FASTQC read quality control
	fastqc: false
	## Run computeGCBias quality control
	gcbias: false
	## Retain only de-duplicated reads/read pairs
	dedup: false
	## Retain only reads with at least the given mapping quality
	mapq: 0
	## Retain only reads mapping in proper pairs
	properpairs: false
	## Mate orientation in paired-end experiments for Bowtie2 mapping
	## (default "--fr" is appropriate for Illumina sequencing)
	mate_orientation: --fr
	## Median/mean fragment length, only relevant for single-end data (default: 200)
	fragment_length: 200
	bowtie_opts:
	qualimap: false
	verbose: False
	################################################################################
	# Call snakemake directly, i.e. without using the wrapper script:
	#
	# Please save a copy of this config yaml file and provide an adjusted config
	# via '--configfile' parameter!
	# example call:
	#
	# snakemake --snakefile /path/to/snakemake_workflows/workflows/DNA-mapping/Snakefile
	#           --configfile /path/to/snakemake_workflows/workflows/DNA-mapping/defaults.yaml
	#           --directory /path/to/outputdir
	#           --cores 32
	################################################################################


Structure of output directory
-------------------------------

The DNA mapping pipeline will generate output of the following structure:

.. code:: bash

    $ tree -d -L 2 output-dir/

.. parsed-literal::

    output-dir/
    |-- Bowtie2
    |-- FASTQ
    |-- FastQC
    |-- Picard_qc
    |   |-- AlignmentSummaryMetrics
    |   |-- InsertSizeMetrics
    |   |-- MarkDuplicates
    |-- Qualimap_qc
    |-- bamCoverage
    |-- cluster_logs
    |-- deepTools_qc
    |   |-- multiBamSummary
    |   |-- plotCorrelation
    |   |-- plotCoverage
    |   `-- plotPCA
    `-- filtered_bam


When enabling read trimming, additional directories will be generated containing the trimmed FASTQ files and, optionally, the FASTQC output on the trimmed FASTQ files.


.. argparse::
   :filename: ../workflows/DNA-mapping/DNA-mapping
   :func: parse_args
   :prog: DNA-mapping
