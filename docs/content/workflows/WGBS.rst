.. _WGBS:

WGBS
============

Input requirements :
-------------------------------------------
This pipeline requires paired-end reads fastq files and a bisulfite converted genome as inputs. 
Optional inputs include bed files with genomic intervals of interest, used to aggregate single CpG values over; a sample sheet with grouping information to use in differential methylation analysis; a blacklist bed file with genomic positions corresponding to known snps to mask single CpG methylation values.


What it does:
-------------------------------------------
Optionally trimmed reads are mapped to reference genome using a bisulfite-specific aligner (bwa-meth).
Quality metrics are collected and synthesized in a QC report, including bisulfite conversion rate, mapping rate, percentage CpGs covered a least 10x, methylation bias.
Methylation ratios are extracted (MethylDackel) for CpG positions in the reference genome and filtered for minimum coverage (10x), snp allelic frequency (<0.25 illegitimate bases).
If sample sheet is provided, logit-transformed beta values for CpG positions are tested for differential methylation using limma.
Metilene is called to detect de novo DMRs. In addition to the nonparametric statistics output by metilene, limma-derived statistics are recalculated for DMRs, which are further annotated with nearest gene information.
If bed file(s) with genomic intervals of interest are provided, methylation ratios are aggregated over those and limma is used on logit-transformed methylation ratios to test for differential methylation.


.. image:: ../images/WGBS_pipeline.png

Configuration file
---------------------

.. code:: bash

    $ cat snakemake_workflows/workflows/WGBS/defaults.yaml

.. parsed-literal::

	################################################################################
	# This file is the default configuration of the WGBS workflow!
	#
	# In order to adjust some parameters, please either use the wrapper script
	# (eg. /path/to/snakemake_workflows/workflows/WGBS/WGBS)
	# or save a copy of this file, modify necessary parameters and then provide
	# this file to the wrapper or snakemake via '--configfile' option
	# (see below how to call the snakefile directly)
	#
	# Own parameters will be loaded during snakefile executiuon as well and hence
	# can be used in new/extended snakemake rules!
	################################################################################
	## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
	wdir:
	configfile:
	cluster_configfile:
	local: False
	max_jobs: 12
	nthreads: 8
	snakemake_options: '--use-conda --conda-prefix /data/processing/conda/envs'
	tempdir: /data/extended/
	## directory with fastq files
	readin:
	## directory with fastqc files for auto quality trimming
	fqcin:
	## Genome information
	genome:
	convrefpath:
	convRef: False
	###list of bed files to process
	intList: []
	###SNP black list (bed file)
	blackList: 
	###sample Info
	sampleInfo:
	###inclusion bounds for methylation extraction
	mbias_ignore: auto
	## FASTQ file extension (default: ".fastq.gz")
	ext: .fastq.gz
	## paired-end read name extension (default: ['_R1', "_R2"])
	reads: [_R1, _R2]
	## Number of reads to downsample from each FASTQ file
	downsample:
	## Options for trimming
	trimReads: auto
	nextera: False
	trimThreshold: 10
	trimOtherArgs:
	verbose: False
	################################################################################
	# Call snakemake directly, i.e. without using the wrapper script:
	#
	# Please save a copy of this config yaml file and provide an adjusted config
	# via '--configfile' parameter!
	# example call:
	#
	# snakemake --snakefile /path/to/snakemake_workflows/workflows/WGBS/Snakefile
	#           --configfile /path/to/snakemake_workflows/workflows/WGBS/defaults.yaml
	#           --directory /path/to/outputdir
	#           --cores 32
	################################################################################


Structure of output directory
--------------------------------

The WGBS pipeline will generate output as follows:

.. code:: bash

    $ tree -d -L 2 output_dir/

.. parsed-literal::

    output_dir
    |-- FASTQ
    |-- FASTQ_Cutadapt
    |   `-- logs
    |-- FASTQ_downsampled
    |   `-- logs
    |-- FastQC_Cutadapt
    |   `-- logs
    |-- FastQC_In
    ...
    |   `-- logs
    |-- QC_metrics
    |   `-- logs
    |-- aggregate_stats_limma
    |   `-- logs
    |-- aux_files
    |   `-- logs
    |-- bams
    |   `-- logs
    |-- cluster_logs
    |-- methXT
    |   `-- logs
    |-- metilene_out
    |   `-- logs
    `-- singleCpG_stats_limma
        `-- logs

Aggregate stats will be calculated if user provides at least one bed file with genomic intervals of interest. Differential methylation analysis or DMR detection will only be run if user provides a sample sheet.

Example output plots 
--------------------------------

Using data from Habibi et al., Cell Stem Cell 2013 corresponding to mouse chr6:4000000-6000000, following plots could be obtained:

.. image:: ../images/limdat.LG.CC.PCA.png

.. image:: ../images/Beta.MeanXgroup.all.violin.png


Argparse
--------------------------------

.. argparse::
   :filename: ../workflows/WGBS/WGBS
   :func: parse_args
   :prog: WGBS
