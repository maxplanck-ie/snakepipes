.. _ChIP-seq:

ChIP-seq
==================

Configuration file
---------------------

.. code:: bash

    $ cat snakemake_workflows/workflows/ChIP-seq/defaults.yaml

.. parsed-literal::

	################################################################################
	# This file is the default configuration of the ChIP-seq workflow!
	#
	# In order to adjust some parameters, please either use the wrapper script
	# (eg. /path/to/snakemake_workflows/workflows/ChIP-seq/ChIP-seq)
	# or save a copy of this file, modify necessary parameters and then provide
	# this file to the wrapper or snakmake via '--configfile' option
	# (see below how to call the snakefile directly)
	#
	# Own parameters will be loaded during snakefile executiuon as well and hence
	# can be used in new/extended snakemake rules!
	################################################################################
	## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
	configfile:
	tempdir: /data/extended/
	local: false
	max_jobs: 5
	snakemake_options:
	## workingdir need to be required DNA-mapping output dir, 'outdir' is set to workingdir internally
	workingdir:
	## preconfigured target genomes (mm9,mm10,dm3,...) , see /path/to/snakemake_workflows/shared/organisms/
	## Value can be also path to your own genome config file!
	genome:
	## paired end data?
	paired: true
	## Bin size of output files in bigWig format
	bw_binsize: 25
	## Median/mean fragment length, only relevant for single-end data (default: 200)
	fragment_length: 200
	verbose: false
	################################################################################
	# Call snakemake directly, i.e. without using the wrapper script:
	#
	# Please save a copy of this config yaml file and provide an adjusted config
	# via '--configfile' parameter!
	# example call:
	#
	# snakemake --snakefile /path/to/snakemake_workflows/workflows/ChIP-seq/Snakefile
	#           --configfile /path/to/snakemake_workflows/workflows/ChIP-seq/defaults.yaml
	#           --directory /path/to/outputdir
	#           --cores 32
	################################################################################

Input requirements
---------------------------

The DNA mapping pipeline generates output that is fully compatible with the ChIP-seq pipeline input requirements!
When running the ChIP-seq pipeline, please specify the output directory of DNA-mapping pipeline as the working directory (-w).

The ChIP-seq pipeline requires at least the following input files for each sample that is specified in the configuration file:

.. code:: bash

    $ tree -L 2 output-dir/filtered_bam/ output-dir/Picard_qc/

..parsed-literal::

    output-dir/filtered_bam/
    |-- my_H3K27ac_sample.filtered.bam
    |-- my_H3K27ac_sample.filtered.bam.bai
    |-- my_H3K27me3_sample.filtered.bam
    |-- my_H3K27me3_sample.filtered.bam.bai
    |-- my_Input_sample.filtered.bam
    `-- my_Input_sample.filtered.bam.bai
    output-dir/Picard_qc/
    |-- AlignmentSummaryMetrics
    |   |-- my_H3K27ac_sample.alignment_summary_metrics.txt
    |   |-- my_H3K27me3_sample.alignment_summary_metrics.txt
    |   `-- my_Input_sample.alignment_summary_metrics.txt
    |-- InsertSizeMetrics
    |   |-- my_H3K27ac_sample.insert_size_metrics.txt
    |   |-- my_H3K27me3_sample.insert_size_metrics.txt
    |   `-- my_Input_sample.insert_size_metrics.txt
    `-- MarkDuplicates
        |-- my_H3K27ac_sample.mark_duplicates_metrics.txt
        |-- my_H3K27me3_sample.mark_duplicates_metrics.txt
        `-- my_Input_sample.mark_duplicates_metrics.txt


Structure of output directory
--------------------------------

The ChIP-seq pipeline will generate additional output as follows:

.. code:: bash

    $ tree -d -L 2 output-dir/

..parsed-literal::

    output-dir/
    ...
    |-- MACS2
    ...
    |-- QC_report
    ...
    |-- cluster_logs
    |-- deepTools_ChIP
    |   |-- bamCompare
    |   |-- plotEnrichment
    |   `-- plotFingerprint
    ...
    `-- histoneHMM

The tool `histoneHMM` will only be run if at least one sample is annotated as broad (IP enrichment).

.. argparse::
   :filename: ../snakePipes/workflows/ChIP-seq/ChIP-seq
   :func: parse_args
   :prog: ChIP-seq
