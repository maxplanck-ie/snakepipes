ChIP-seq workflow v0.3.2.1 - MPI-IE workflow for ChIP-seq analysis
Fabian Kilpert, Andreas Richter
June 27, 2016


## Quick start ##

A **typical ChIP-seq analysis of human samples** starting from paired-end FASTQ files in the directory `input-dir`:

    $ ls /path/to/input-dir/
    my_H3K27ac_sample_R1.fastq.gz  my_H3K27me3_sample_R1.fastq.gz  my_Input_sample_R1.fastq.gz
    my_H3K27ac_sample_R2.fastq.gz  my_H3K27me3_sample_R2.fastq.gz  my_Input_sample_R2.fastq.gz
    $
    $ /data/manke/repository/scripts/snakemake_workflows/snakemake_workflows_0.3.2.1/DNA-mapping hs37d5 \
          -i /path/to/input-dir -o /path/to/output-dir --dedup && \
      /data/manke/repository/scripts/snakemake_workflows/snakemake_workflows_0.3.2.1/ChIP-seq hs37d5 chip-seq.config.yaml \
          -d /path/to/outputdir

All individual jobs of the workflow will be submitted to the Slurm queue. To run the workflow locally, use the parameter `--local` for local mode and the parameter `-j 48` to specify the number of available CPU cores (here: 48).

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

    $ /data/manke/repository/scripts/snakemake_workflows/snakemake_workflows_0.3.2.1/DNA-mapping -h


### Configuration file ###

    $ cat /data/manke/repository/scripts/snakemake_workflows/snakemake_workflows_0.3.2.1/workflows/DNA-mapping/example.config.yaml

    ################################################################################
    # Please comment out and adjust the following configuration variables if their
    # default values are not appropriate for your input data
    ################################################################################
    ## paired-end read name extension (default: ["_R1", "_R2"])
    # reads: ["_R1", "_R2"]

    ## FASTQ file extension (default: ".fastq.gz")
    # ext: ".fastq.gz"

    ## Mate orientation in paired-end experiments for Bowtie2 mapping
    ## (default "--fr" is appropriate for Illumina sequencing)
    # mate_orientation: "--fr"

    ## Median/mean fragment length, only relevant for single-end data (default: 200)
    # fragment_length: 200

    ## Options to TrimGalore (default: "--stringency 2")
    # trim_galore_opts: "--stringency 2"

    ################################################################################
    # The following configuration variables can only be set via this configuration
    # file when calling snakemake directly, i.e. without using the wrapper script
    # DNA-mapping. When calling the script DNA-mapping, the values set below are
    # ignored.
    # example call: snakemake --snakefile /path/to/snakemake_workflows/workflows/DNA-mapping/Snakefile
    #                         --cores 32
    #                         --configfile /path/to/snakemake_workflows/workflows/DNA-mapping/example.config.yaml
    #                         --directory /path/to/outputdir
    ################################################################################
    ## Snakemake workflow directory
    # maindir: /path/to/snakemake_workflows/
    ## Input directory with FASTQ files
    # indir: /path/to/inputdir/
    ## Genome acronym
    # genome: hs37d5
    ## Number of reads to downsample from the head of each FASTQ file
    # downsample: 50000
    ## Trim reads with TrimGalore
    # trim: True
    ## Run FASTQC read quality control
    # fastqc: True
    ## Retain only de-duplicated reads/read pairs
    # dedup: True
    ## Retain only reads mapping in proper pairs
    # properpairs: True
    ## Retain only reads with at least the given mapping quality
    # mapq: 1
    ## Bin size of output files in bigWig format
    # bw_binsize: 10
    ## Run computeGCBias quality control
    # gcbias: False
    ##  Verbose output
    # verbose: True


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

When enabling read trimming with TrimGalore, additional directories will be generated containing the trimmed FASTQ files and, optionally, the FASTQC output on the trimmed FASTQ files.


## ChIP-seq pipeline ##
### Usage ###

    $ /data/manke/repository/scripts/snakemake_workflows/snakemake_workflows_0.3.2.1/ChIP-seq -h


### Configuration file ###

    $ cat /data/manke/repository/scripts/snakemake_workflows/snakemake_workflows_0.3.2.1/workflows/ChIP-seq/example.config.yaml
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

    ################################################################################
    # Please comment out and adjust the following configuration variables if their
    # default values are not appropriate for your input data
    # ################################################################################
    ## Median/mean fragment length, only relevant for single-end data (default: 200)
    # fragment_length: 200

    ################################################################################
    # The following configuration variables can only be set via this configuration
    # file when calling snakemake directly, i.e. without using the wrapper script
    # ChIP-seq. When calling the script ChIP-seq, the values set below are ignored.
    # example call: snakemake --snakefile /path/to/snakemake_workflows/workflows/ChIP-seq/Snakefile
    #                         --cores 32
    #                         --configfile /path/to/snakemake_workflows/workflows/ChIP-seq/example.config.yaml
    #                         --directory /path/to/outputdir
    ################################################################################
    ## Snakemake workflow directory
    # maindir: /path/to/snakemake_workflows/
    ## Genome acronym
    # genome: hs37d5
    ## Input data is paired-end
    # paired: True
    ## Bin size of output files in bigWig format
    # bw_binsize: 10
    ##  Verbose output
    # verbose: True


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

Further organisms can be supported by adding a genome configuration file `my_organism.py` in the following style to the `snakemake_workflows` directory:

    $ cat /data/manke/repository/scripts/snakemake_workflows/snakemake_workflows_0.3.2.1/shared/organisms/hs37d5.py

    genome_size = 2900338458
    genome_fasta = "/data/repository/organisms/hs37d5_ensembl/genome_fasta/genome.fa"
    genome_index = "/data/repository/organisms/hs37d5_ensembl/genome_fasta/genome.fa.fai"
    genome_2bit = "/data/repository/organisms/hs37d5_ensembl/genome_fasta/genome.2bit"
    bowtie2_index = "/data/repository/organisms/hs37d5_ensembl/BowtieIndex/genome"
    genes_bed = "/data/repository/organisms/hs37d5_ensembl/gencode/release_19/genes.bed"
    genes_gtf = "/data/repository/organisms/hs37d5_ensembl/gencode/release_19/genes.gtf"
    blacklist_bed = "/data/repository/organisms/hs37d5_ensembl/ENCODE/hs37d5_extended_Encode-blacklist.bed"

If no blacklist regions are available for your organism of interest, please set `blacklist_bed = False`.
