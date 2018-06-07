.. _WGBS:

WGBS
============

Input requirements and outputs:
-------------------------------------------
This pipeline requires paired-end reads fastq files and a bisulfite converted genome as inputs. Methylation values per CpG dinucleotide are extracted after read mapping. 
Statistical analysis of differential methylation is performed if user provides a bed file with genomic intervals to aggregate methylation values over. De novo DMR calling, statistical evaluation and annotation with nearest gene are performed automatically. Various quality metrics are collected and a QC report is output.

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

The WGBS pipeline will generate additional output as follows:

.. code:: bash

    $ tree -d -L 2 output-dir/

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


Usage example
--------------------------------

.. parsed-literal::

python3 snakemake_workflows/workflows/WGBS/WGBS --indir /data/manke/group/snakepipes_test_files/WGBS/reads --wdir /data/processing3/WGBS_snakepipe_example_OUT_logging_auto --intList /data/manke/group/snakepipes_test_files/WGBS/danRer10.cpgIsland.ext.sorted.chr25.noCHR.bed --sampleInfo /data/manke/group/snakepipes_test_files/WGBS/example_sampleSheet.csv GRCz10


Argparse
--------------------------------

.. argparse::
   :filename: ../workflows/WGBS/WGBS
   :func: parse_args
   :prog: WGBS
