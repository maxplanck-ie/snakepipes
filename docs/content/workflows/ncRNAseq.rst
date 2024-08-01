.. _ncRNAseq:

ncRNAseq
=================

What it does
------------

The snakePipes ncRNAseq workflow allows users to process their single or paired-end
ribosomal-depleted RNA-seq fastq files upto the point of gene/transcript/repeat-element counts and differential expression.
Repeat elements are quantified and tested for differential expression at the name, family and class level. Since changes in repeat element expression tend to be unidirectional, size factors from gene expression are used when normalizing repeat element expression.

Note that in addition to the normal GTF file, this pipeline requires a repeat masker output file, which can be downloaded from UCSC or other sites. The chromosome names here must match that used in the other files.

.. image:: ../images/ncRNAseq_pipeline.png

Input requirements
------------------

The only requirement is a directory of gzipped fastq files. Files could be single or paired end, and the read extensions could be modified using the keys in the ``defaults.yaml`` file below.

.. _ncRNAconfig:

Configuration file
~~~~~~~~~~~~~~~~~~

There is a configuration file in ``snakePipes/workflows/ncRNAseq/defaults.yaml``::

    ## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
    pipeline: ncRNAseq
    outdir:
    configFile:
    clusterConfigFile:
    local: False
    maxJobs: 5
    ## directory with fastq files
    indir:
    ## preconfigured target genomes (mm9,mm10,dm3,...) , see /path/to/snakemake_workflows/shared/organisms/
    ## Value can be also path to your own genome config file!
    genome:
    ## FASTQ file extension (default: ".fastq.gz")
    ext: '.fastq.gz'
    ## paired-end read name extension (default: ["_R1", "_R2"])
    reads: ["_R1","_R2"]
    ## assume paired end reads
    pairedEnd: True
    ## Number of reads to downsample from each FASTQ file
    downsample:
    ## Options for trimming
    trim: False
    trimmer: fastp
    trimmerOptions:
    ## further options
    mode: alignment,deepTools_qc
    sampleSheet:
    bwBinSize: 25
    fastqc: False
    fragmentLength: 200
    libraryType: 2
    ## supported mappers: STAR HISAT2
    aligner: STAR
    alignerOptions: "--outSAMstrandField intronMotif --outFilterMultimapNmax 1000 --outFilterMismatchNoverLmax 0.1 --outSAMattributes Standard"
    verbose: False
    plotFormat: png
    #### Flag to control the pipeline entry point
    fromBAM: False
    bamExt: '.bam'
    #umi_tools
    UMIBarcode: False
    bcPattern: NNNNCCCCCCCCC #default: 4 base umi barcode, 9 base cell barcode (eg. RELACS barcode)
    UMIDedup: False
    UMIDedupSep: "_"
    UMIDedupOpts: --paired


Apart from the common workflow options (see :ref:`running_snakePipes`), the following parameters are useful to consider:

* **aligner**: The only choice at the moment is `STAR <https://github.com/alexdobin/STAR>`__. 

* **alignerOptions**: Options to pass on to your chosen aligner. Note that library type and junction definitions don't have to be passed to the aligners as options, as they are handeled either automatically, or via other parameters.

* **plotFormat**: You can switch the type of plot produced by all deeptools modules using this option. Possible choices : png, pdf, svg, eps, plotly


Differential expression
-----------------------

Like the other workflows, differential expression can be performed using the ``--sampleSheet`` option and supplying a sample sheet like that below::

    name    condition
    sample1      eworo
    sample2      eworo
    SRR7013047      eworo
    SRR7013048      OreR
    SRR7013049      OreR
    SRR7013050      OreR

.. note:: The first entry defines which group of samples are control. This way, the order of comparison and likewise the sign of values can be changed. The DE analysis might fail if your sample names begin with a number. So watch out for that!

Complex designs with blocking factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the user provides additional columns between 'name' and 'condition' in the sample sheet, the variables stored there will be used as blocking factors in the order they appear in the sample sheet. Eg. if the first line of your sample sheet looks like 'name	batch	condition', this will translate into a formula ``batch + condition``. 'condition' has to be the final column and it will be used for any statistical inference.

Multiple pairwise comparisons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The user may specify multiple groups of independent comparisons by providing a 'group' column after the 'condition' column. This will cause the sample sheet to be split into the groups defined in this column, and a corresponding number of independent pairwise comparisons will be run, one for each split sheet, and placed in separate output folders named accordingly. This will be applied to DESeq2 pairwise comparison.
Specifying a value of 'All' in the 'group' column will cause that sample group to be used in all pairwise comparisons, e.g. if the same set of controls should be used for several different treatment groups.

An example sample sheet with the group information provided looks like this:

name	condition	group
sample1	Control		All
sample2	Control		All
sample3	Treatment	Group1
sample4	Treatment	Group1
sample5	Treatment	Group2
sample6	Treatment	Group2

Analysis modes
--------------

Following analysis (**modes**) are possible using the ncRNAseq workflow:

"alignment"
~~~~~~~~~~~

In this mode,
the pipeline uses **STAR** to create BAM files and **TEtranscripts** to quantify genes and repeat elements.
Differential expression of genes and repeat elements is then performed with **DESeq2**.

"deepTools_qc"
~~~~~~~~~~~~~~

The pipeline provides multiple quality controls through deepTools, which can be triggered
using the **deepTools_qc** mode. It's a very useful add-on with any of the other modes.

.. note:: Since most deeptools functions require an aligned (BAM) file, the deepTools_qc mode will additionally perform the alignment of the fastq files. However this would not interfere with operations of the other modes.

Understanding the outputs
-------------------------

Assuming the pipline was run with ``--mode 'alignment,deepTools_qc'`` on a set of FASTQ files, the structure of the output directory would look like this (files are shown only for one sample) ::

    ├── bamCoverage
    │   ├── sample1.coverage.bw
    │   ├── sample1.RPKM.bw
    │   ├── sample1.uniqueMappings.fwd.bw
    │   ├── sample1.uniqueMappings.rev.bw
    ├── cluster_logs
    ├── deepTools_qc
    │   ├── bamPEFragmentSize
    │   │   ├── fragmentSize.metric.tsv
    │   │   └── fragmentSizes.png
    │   ├── estimateReadFiltering
    │   │   └── sample1_filtering_estimation.txt
    │   ├── logs
    │   │   ├── bamPEFragmentSize.err
    │   │   ├── bamPEFragmentSize.out
    │   │   ├── multiBigwigSummary.err
    │   │   └── plotCorrelation_pearson.err
    │   ├── multiBigwigSummary
    │   │   └── coverage.bed.npz
    │   ├── plotCorrelation
    │   │   ├── correlation.pearson.bed_coverage.heatmap.png
    │   │   ├── correlation.pearson.bed_coverage.tsv
    │   │   ├── correlation.spearman.bed_coverage.heatmap.png
    │   │   └── correlation.spearman.bed_coverage.tsv
    │   ├── plotEnrichment
    │   │   ├── plotEnrichment.png
    │   │   └── plotEnrichment.tsv
    │   └── plotPCA
    │       ├── PCA.bed_coverage.png
    │       └── PCA.bed_coverage.tsv
    ├── DESeq2_sampleSheet
    │   ├── DESeq2_report_genes.html
    │   ├── DESeq2_report_repeat_class.html
    │   ├── DESeq2_report_repeat_family.html
    │   ├── DESeq2_report_repeat_name.html
    │   ├── DESeq2.session_info.txt
    │   ├── genes_counts_DESeq2.normalized.tsv
    │   ├── genes_DEresults.tsv
    │   ├── genes_DESeq.Rdata
    │   ├── repeat_class_counts_DESeq2.normalized.tsv
    │   ├── repeat_class_DEresults.tsv
    │   ├── repeat_class_DESeq.Rdata
    │   ├── repeat_family_counts_DESeq2.normalized.tsv
    │   ├── repeat_family_DEresults.tsv
    │   ├── repeat_family_DESeq.Rdata
    │   ├── repeat_name_counts_DESeq2.normalized.tsv
    │   ├── repeat_name_DEresults.tsv
    │   └── repeat_name_DESeq.Rdata
    ├── FASTQ
    ├── filtered_bam
    │   ├── sample1.filtered.bam
    │   ├── sample1.filtered.bam.bai
    ├── multiQC
    ├── STAR
    └── TEcount
        └── sample1.cntTable


.. note:: The ``_sampleSheet`` suffix for the ``DESeq2_sampleSheet`` is drawn from the name of the sample sheet you use. So if you instead named the sample sheet ``mySampleSheet.txt`` then the folder would be named ``DESeq2_mySampleSheet``. This facilitates using multiple sample sheets.

Apart from the common module outputs (see :ref:`running_snakePipes`), the workflow would produce the following folders:

* **bamCoverage**: This would contain the bigWigs produced by deepTools `bamCoverage <https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html>`__ . Files with suffix ``.coverage.bw`` are raw coverage files, while the files with suffix ``RPKM.bw`` are `RPKM-normalized <https://www.biostars.org/p/273537/>`__ coverage files.

* **deepTools_QC**: (produced in the mode *deepTools_QC*) This contains the quality checks specific for mRNAseq, performed via deepTools. The output folders are names after various deepTools functions and the outputs are explained under `deepTools documentation <deeptools.readthedocs.io>`__. In short, they show the insert size distribution(**bamPEFragmentSize**), mapping statistics (**estimateReadFiltering**), sample-to-sample correlations and PCA (**multiBigwigSummary, plotCorrelation, plotPCA**), and read enrichment on various genic features (**plotEnrichment**)

* **DESeq2_[sampleSheet]**: (produced only if a sample-sheet is provided.) The folder contains the HTML result reports **DESeq2_report_genes.html**, **DESeq2_report_repeat_name.html**, **DESeq2_report_repeat_class.html** and **DESeq2_report_repeat_family.html** as we as the annotated output file from DESeq2 (**genes_DEresults.tsv**, etc.) and normalized counts for all samples, produced via DEseq2 (**genes_counts_DESeq2.normalized.tsv**, etc.) as well as an Rdata file (**genes_DESeq.Rdata**, etc.) with the R objects ``dds <- DESeq2::DESeq(dds)`` and ``ddr <- DDESeq2::results(dds,alpha = fdr)``. Sample name to plotting shape mapping on the PCA plot is limited to 36 samples and skipped otherwise.

* **filtered_bam**: This contains sorted and indexed BAM files that have been filtered to remove secondary alignments. This are used by deepTools and are appropriate for use in IGV.

* **multiQC**: This folder contains the report produced by MultiQC and summarizes alignment metrics from STAR and possibly various deepTools outputs.

* **STAR**: This would contain the output logs of RNA-alignment by STAR. The actual BAM files are removed at the end of the pipeline as they're not compatible with typical visualization programs

* **TEcount**: (produced in the **alignment** mode) This contains the counts files and logs from the TEcount program in the TEtranscripts package. These are used by DESeq2 for differential expression.


Command line options
--------------------

.. argparse::
    :func: parse_args
    :filename: ../snakePipes/workflows/ncRNAseq/ncRNAseq
    :prog: ncRNAseq
    :nodefault:
