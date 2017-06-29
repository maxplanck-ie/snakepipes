MPI-IE Snakemake workflows

Fabian Kilpert, Andreas Richter, Steffen Heyne

2016-2017

## Workflows ##

- DNA-mapping
- ChIP-seq
- RNA-seq
- scRNA-seq

## Quick start ##

A **typical ChIP-seq analysis of human samples** starting from paired-end FASTQ files in the directory `input-dir`:

    $ ls /path/to/input-dir/
    my_H3K27ac_sample_R1.fastq.gz  my_H3K27me3_sample_R1.fastq.gz  my_Input_sample_R1.fastq.gz
    my_H3K27ac_sample_R2.fastq.gz  my_H3K27me3_sample_R2.fastq.gz  my_Input_sample_R2.fastq.gz
    $
    $ snakemake_workflows/DNA-mapping hs37d5 \
          -i /path/to/input-dir -o /path/to/output-dir --dedup && \
      snakemake_workflows/ChIP-seq hs37d5 chip-seq.config.yaml \
          -d /path/to/outputdir

All individual jobs of the workflow will be submitted to the Slurm queue. To run the workflow locally, use the parameter `--local` for local mode and the parameter `-j 48` to specify the maximal number of used CPU threads (here: 48) or concurrent running Slurm jobs (actual used threads are defined in each rule).

A **configuration file is required for the ChIP-seq workflow** and should adhere to the following style (IMPORTANT: Use only whitespace, but NO TAB indentation in this file:

    $ cat chip-seq.config.yaml
    ################################################################################
    # Please specify all ChIP samples plus their matching control/chromatin input
    # sample.
    # Specify for each ChIP sample whether the IP target results in broad/mixed-type
    # enrichment (most histone marks, e.g. H3K4me1, H3K36me3, H3K9me3, H3K27me3)
    # or not. In the latter case, the enrichment is instead punctuate/narrow
    # (e.g. TFs, active histone marks as H3K27ac or H3K4me3).
    #
    # IMPORTANT: Use only whitespace, but NO TAB indentation in this YAML file!
    ################################################################################
    chip_dict:
      my_H3K27ac_sample:
        control: my_Input_sample
        broad: False
      my_H3K27me3_sample:
        control: my_Input_sample
        broad: True


## DNA mapping pipeline ##
### Usage ###

    $ snakemake_workflows/DNA-mapping -h


### Configuration file ###

    $ cat snakemake_workflows/workflows/DNA-mapping/defaults.yaml

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


### Structure of output directory ###

The DNA mapping pipeline will generate output of the following structure:

    $ tree -d -L 2 output-dir/

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


## ChIP-seq pipeline ##
### Usage ###

    $ snakemake_workflows/ChIP-seq -h


### Configuration file ###

    $ cat snakemake_workflows/workflows/ChIP-seq/defaults.yaml

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

### Input requirements ###

The DNA mapping pipeline generates output that is fully compatible with the ChIP-seq pipeline input requirements!

The ChIP-seq pipeline requires at least the following input files for each sample that is specified in the configuration file:

    $ tree -L 2 output-dir/filtered_bam/ output-dir/Picard_qc/

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


### Structure of output directory ###

The ChIP-seq pipeline will generate additional output as follows:

    $ tree -d -L 2 output-dir/

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


## Genome configuration file ##

Further organisms can be supported by adding a genome configuration file `my_organism.yaml` in the following style to the `snakemake_workflows` directory:

    $ cat snakemake_workflows/shared/organisms/hs37d5.yaml

	genome_size: 2900338458
	genome_fasta: "/SOMEPATH/hs37d5_ensembl/genome_fasta/genome.fa"
	genome_index: "/SOMEPATH/hs37d5_ensembl/genome_fasta/genome.fa.fai"
	genome_2bit: "/SOMEPATH/hs37d5_ensembl/genome_fasta/genome.2bit"
	bowtie2_index: "/SOMEPATH/hs37d5_ensembl/BowtieIndex/genome"
	hisat2_index: "/SOMEPATH/hs37d5_ensembl/HISAT2Index/genome"
	known_splicesites: "/SOMEPATH/hs37d5_ensembl/gencode/release_19/HISAT2/splice_sites.txt"
	star_index: "/SOMEPATH/hs37d5_ensembl/STARIndex/"
	genes_bed: "/SOMEPATH/hs37d5_ensembl/gencode/release_19/genes.bed"
	genes_gtf: "/SOMEPATH/hs37d5_ensembl/gencode/release_19/genes.gtf"
	blacklist_bed: "/SOMEPATH/hs37d5_ensembl/ENCODE/hs37d5_extended_Encode-blacklist.bed"

If no blacklist regions are available for your organism of interest, leave it empty `blacklist_bed: `
