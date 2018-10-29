.. _DNA-mapping:

DNA-mapping
===========

What it does
------------

This is the primary DNA-mapping pipeline. It can be used both alone or upstream of the ATAC-seq and ChIP-seq pipelines. This has a wide array of options, including trimming and various QC steps (e.g., marking duplicates and plotting coverage and PCAs). In addition, basic coverage tracks are created to facilitate viewing the data in IGV.

.. image:: ../images/DNAmapping_pipeline.png

Input requirements
------------------

The only requirement is a directory of gzipped fastq files. Files could be single or paired end, and the read extensions could be modified using the keys in the `defaults.yaml` file below.

Configuration file
~~~~~~~~~~~~~~~~~~

There is a configuration file in `snakePipes/workflows/DNA-mapping/defaults.yaml`::

    ## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
    pipeline: dna-mapping
    outdir:
    configfile:
    cluster_configfile:
    local: False
    max_jobs: 5
    ## directory with fastq files
    indir:
    ## preconfigured target genomes (mm9,mm10,dm3,...) , see /path/to/snakemake_workflows/shared/organisms/
    ## Value can be also path to your own genome config file!
    genome:
    ## FASTQ file extension (default: ".fastq.gz")
    ext: '.fastq.gz'
    ## paired-end read name extension (default: ['_R1', "_R2"])
    reads: [_R1, _R2]
    ## mapping mode
    mode: mapping
    mapping_prg: Bowtie2
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
    ## other Bowtie2 stuff
    insert_size_max: 1000
    bowtie_opts:
    plot_format: png
    ## Median/mean fragment length, only relevant for single-end data (default: 200)
    fragment_length: 200
    qualimap: false
    verbose: false

Many of these options can be more conveniently set on the command-line (e.g., `--qualimap` sets `qualimap: true`). However, you may need to change the `reads:` setting if your paired-end files are not denoted by `sample_R1.fastq.gz` and `sample_R2.fastq.gz`, but rather `sample_1.fastq.gz` and `sample_2.fastq.gz`.

Understanding the outputs
--------------------------

The DNA mapping pipeline will generate output of the following structure::

    .
    ├── bamCoverage
    ├── Bowtie2
    ├── deepTools_qc
    │   ├── bamPEFragmentSize
    │   ├── estimateReadFiltering
    │   ├── multiBamSummary
    │   ├── plotCorrelation
    │   ├── plotCoverage
    │   └── plotPCA
    ├── FASTQ
    ├── FastQC
    ├── filtered_bam
    ├── multiQC
    │   └── multiqc_data
    ├── Picard_qc
    │   └── AlignmentSummaryMetrics
    └── Sambamba

In addition to the FASTQ module results (see :doc:`running_snakePipes`), the workflow produces the following outputs:

 * **Bowtie2** : Contains the BAM files after mapping with `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`__ and indexed by `Samtools <http://www.htslib.org/>`__.

 * **filtered_bam** : Contains the BAM files filtered by the provided criteria, such as mapping quality (`--mapq`) or PCR duplicates (`--dedup`). This file is used for most downstream analysis in the DNA-mapping and ChIP-seq/ATAC-seq pipeline.

 * **bamCoverage** : Contains the coverage files (`bigWig format <https://genome.ucsc.edu/goldenpath/help/bigWig.html>`__) produced from the BAM files by `deepTools bamCoverage <https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html>`__ . The files are either raw, or 1x normalized (by sequencing depth). They are useful for plotting and inspecting the data in IGV.

 * **deepTools_qc** : Contains various QC files and plots produced by deepTools on the filtered BAM files. These are very useful for evaluation of data quality. The folders are named after the tools. Please look at the `deepTools documentation <https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html>`__ on how to interpret the outputs from each tool.

* **Picard_qc** : Contains the output by `picard CollectAlignmentSummaryMetrics <https://broadinstitute.github.io/picard/command-line-overview.html>`__ tool. This output is used for the evaluation of reads within peaks by the ChIp-Seq and ATAC-seq workflows.

* **Sambamba** : Contains the alignment metrices evaluated on the BAM files by `Sambamba <http://lomereiter.github.io/sambamba/>`__.

A number of other directories may optionally be present if you specified read trimming, using Qualimap, or a variety of other options. These are typically self-explanatory.

A fair number of useful QC plots are or can be generated by the pipeline. These include correlation and PCA plots as well as the output from MultiQC.

.. image:: ../images/DNAmapping_correlation.png

Command line options
--------------------

.. argparse::
   :func: parse_args
   :filename: ../snakePipes/workflows/DNA-mapping/DNA-mapping
   :prog: DNA-mapping
   :nodefault:
