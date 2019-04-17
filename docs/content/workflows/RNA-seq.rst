.. _RNA-seq:

RNA-seq
=======

What it does
------------

The snakePipes RNA-seq workflow allows users to process their single or paired-end
RNA-Seq fastq files upto the point of gene/transcript-counts and differential expression.
It also allows full allele-specific RNA-Seq analysis (up to allele-specific
differential expression) using the *allelic-mapping* mode.

.. image:: ../images/RNAseq_pipeline.png

Input requirements
------------------

The only requirement is a directory of gzipped fastq files. Files could be single or paired end, and the read extensions could be modified using the keys in the ``defaults.yaml`` file below.

.. _RNAconfig:

Configuration file
~~~~~~~~~~~~~~~~~~

There is a configuration file in ``snakePipes/workflows/RNA-seq/defaults.yaml``::

    ## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
    pipeline: rna-seq
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
    ## paired-end read name extension (default: ["_R1", "_R2"])
    reads: ["_R1","_R2"]
    ## Number of reads to downsample from each FASTQ file
    downsample:
    ## Options for trimming
    trim: False
    trim_prg: cutadapt
    trim_options:
    ## further options
    mode: alignment-free,deepTools_qc
    sampleSheet:
    bw_binsize: 25
    fastqc: False
    featurecounts_options: -C -Q 10 --primary
    filter_annotation:
    fragment_length: 200
    library_type: 2
    salmon_index_options: --type quasi -k 31
    dnaContam: False
    ## supported mappers: STAR HISAT2
    mapping_prg: STAR
    ## N.B., setting --outBAMsortingBinsN too high can result in cryptic errors
    star_options: --outBAMsortingBinsN 30
    hisat_options:
    verbose: False
    plot_format: png
    # for allele-spcific mapping
    snp_file:
    Nmasked_index:


Apart from the common workflow options (see :ref:`running_snakePipes`), the following parameters are useful to consider:

* **mapping_prg**: You can choose either `STAR <https://github.com/alexdobin/STAR>`__ or `HISAT2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`__. While HISAT2 mapping is usually faster than STAR, we keep STAR as the default aliger due to it's superior accuracy (see `this paper <https://www.nature.com/nmeth/journal/v14/n2/abs/nmeth.4106.html>`__).

* **star_options/hisat_options**: Options to pass on to your chosen aligner. Note that library type and junction definitions doesn't have to be passed to the aligners as options, as they are handeled either automatically, or via other parameters.

* **featurecounts_options**: Options to pass to featureCounts (in case the ``alignment`` or ``allelic-mapping`` mode is used). Note that the paired-end information is automatically passed to featurecounts via the workflow, and the summerization is always performed at **gene level**, since the workflow assumes that featurecounts output is being used for gene-level differential expression analysis. If you wish to perform a **transcript-level** DE analysis, please choose the mode **alignment-free**.

* **filter_annotation**: Options you can pass on to filter the original `GTF file <http://genome.ucsc.edu/goldenPath/help/customTrack.html#GTF>`__. This is useful in case you want to filter certain kind of transcripts (such as pseudogenes) before running the counts/DE analysis.

* **library_type**: The default library-type is suitable for most RNAseq protocols (using Illumina `Tru-Seq <https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-rna-v2.html>`__). Change this option in case you have a different strandednes.

* **salmon_index_options**: In the ``alignment-free`` mode (see below), this option allows you to change the type of index created by salmon. New users can leave it to default.

* **dnaContam**: Enable this to test for possible DNA contamination in your RNA-seq samples. DNA contamination is quantified as the fraction of reads falling into intronic and intergenic regions, compared to those falling into exons. Enabling this option would produce a directory called ``GenomicContamination`` with ``.tsv`` files containing this information.

* **plot_format**: You can switch the type of plot produced by all deeptools modules using this option. Possible choices : png, pdf, svg, eps, plotly

* **snp_file**: For the ``allelic-mapping`` mode. The ``snp_file`` is the file produced by `SNPsplit <https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/>`__ after creating a dual-hybrid genome. The file has the suffix ``.vcf``.

* **Nmasked_index**: For the ``allelic-mapping`` mode. The ``Nmasked_index`` refers to the **basename** of the index file created using STAR, on the SNPsplit output.

.. note:: snp_file and Nmasked_index file could be specified in case you already have this and would like to re-run the analysis on new data. In case you are running the allele-specific analysis for the first time, you would need a VCF file and the name of the two strains. In this case the ``snp_file`` as well as the ``Nmasked_index`` files would be automatically created by the workflow using SNPsplit.


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

If the user provides additional columns between 'name' and 'condition' in the sample sheet, the variables stored there will be used as blocking factors in the order they appear in the sample sheet. Condition will be the final column and it will be used for any statistical inference. 

Analysis modes
--------------

Following analysis (**modes**) are possible using the RNA-seq workflow:

"alignment"
~~~~~~~~~~~

In this mode,
the pipeline uses one of the selected aligners to create BAM files, followed by
gene-level quantification using **featureCounts**. Gene-level differential expression
analysis is then performed using **DESeq2**.

"allelic-mapping"
~~~~~~~~~~~~~~~~~

**allelic-mapping** mode follows a similar process as the "mapping" mode, however the
alignment performed on an allele-masked genome, followed by allele-specific splitting
of mapped files. Gene-level quantification is performed for each allele using **featureCounts**.
Allele-specific, gene-level differential expression analysis is then performed using **DESeq2**.

.. note:: **allelic-mapping** mode is mutually exclusive with **mapping** mode

"alignment-free"
~~~~~~~~~~~~~~~~

In this mode,
the pipeline uses **salmon** to perform transcript-level expression quantification. This mode
performs both transcript-level differential expression (using **Sleuth**), and gene-level
differential expression (using **wasabi**, followed by **DESeq2**).

.. note:: The salmon index is recreated each time in alignment-free mode. This is done to facilitate changing how the GTF file is filtered, which necessitates reindexing.

"deepTools_qc"
~~~~~~~~~~~~~~

The pipeline provides multiple quality controls through deepTools, which can be triggered
using the **deepTools_qc** mode. It's a very useful add-on with any of the other modes.

.. note:: Since most deeptools functions require an aligned (BAM) file, the deepTools_qc mode will additionally perform the alignment of the fastq files. However this would not interfere with operations of the other modes.

Understanding the outputs
---------------------------

Assuming the pipline was run with ``--mode 'alignment-free,alignment,deepTools_qc'`` on a set of FASTQ files, the structure of the output directory would look like this (files are shown only for one sample) ::

    ├── Annotation
    │   ├── filter_command.txt
    │   ├── genes.annotated.bed
    │   ├── genes.filtered.bed
    │   ├── genes.filtered.fa
    │   ├── genes.filtered.gtf
    │   ├── genes.filtered.symbol
    │   ├── genes.filtered.t2g
    ├── bamCoverage
    │   ├── logs
    │   ├── sample1.coverage.bw
    │   ├── sample1.RPKM.bw
    │   ├── sample1.uniqueMappings.fwd.bw
    │   └── sample1.uniqueMappings.rev.bw
    ├── cluster_logs
    ├── deepTools_qc
    │   ├── bamPEFragmentSize
    │   │   ├── fragmentSize.metric.tsv
    │   │   └── fragmentSizes.png
    │   ├── estimateReadFiltering
    │   │   └── sample1_filtering_estimation.txt
    │   ├── logs
    │   ├── multiBigwigSummary
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
    ├── DESeq2_Salmon_sampleSheet
    │   ├── DESeq2_Salmon.err
    │   ├── DESeq2_Salmon.out
	│   ├── citations.bib
	│   ├── DESeq2_report_files
	│   ├── DESeq2_report.html
	│   ├── DESeq2_report.Rmd
	│   ├── DESeq2.session_info.txt
	│   ├── DEseq_basic_counts_DESeq2.normalized.tsv
	│   ├── DEseq_basic_DEresults.tsv
	│   └── DEseq_basic_DESeq.Rdata
    ├── DESeq2_sampleSheet
    │   ├── DESeq2.err
    │   ├── DESeq2.out
	│   ├── citations.bib
	│   ├── DESeq2_report_files
	│   ├── DESeq2_report.html
	│   ├── DESeq2_report.Rmd
	│   ├── DESeq2.session_info.txt
	│   ├── DEseq_basic_counts_DESeq2.normalized.tsv
	│   ├── DEseq_basic_DEresults.tsv
	│   └── DEseq_basic_DESeq.Rdata
    ├── FASTQ
    │   ├── sample1_R1.fastq.gz
    │   └── sample1_R2.fastq.gz
	├── featureCounts
	│   ├── counts.tsv
	│   ├── sample1.counts.txt
	│   ├── sample1.counts.txt.summary
	│   ├── sample1.err
	│   ├── sample1.out
    ├── multiQC
    │   ├── multiqc_data
    │   ├── multiQC.err
    │   ├── multiQC.out
    │   └── multiqc_report.html
	├── QC_report
	│   └── QC_report_all.tsv
    ├── RNA-seq.cluster_config.yaml
    ├── RNA-seq.config.yaml
	├── RNA-seq_organism.yaml
    ├── RNA-seq_pipeline.pdf
    ├── RNA-seq_run-1.log
    ├── Salmon
    │   ├── counts.genes.tsv
    │   ├── counts.tsv
    │   ├── Salmon_counts.log
    │   ├── Salmon_genes_counts.log
    │   ├── Salmon_genes_TPM.log
    │   ├── SalmonIndex
    │   ├── Salmon_TPM.log
    │   ├── sample1
    │   ├── sample1.quant.genes.sf
    │   ├── sample1.quant.sf
    │   ├── TPM.genes.tsv
    │   └── TPM.tsv
	├── sleuth_Salmon_sampleSheet
	│   ├── logs
	│   ├── MA-plot.pdf
	│   ├── sleuth_live.R
	│   ├── so.rds
	│   └── Wald-test.results.tsv
    └── STAR
  		├── logs
        ├── sample1
        ├── sample1.bam
        └── sample1.bam.bai

Apart from the common module outputs (see :ref:`running_snakePipes`), the workflow would produce the following folders:

* **Annotation**: This folder would contain the GTF and BED files used for analysis. In case the file has been filtered using the ``--filter_annotation`` option (see :ref:`RNAconfig`), this would contain the filtered files.

* **STAR/HISAT2**: (not produced in mode *alignment-free*) This would contain the output of RNA-alignment by STAR or HISAT2 (indexed `BAM files <http://samtools.github.io/hts-specs/SAMv1.pdf>`__).

* **featureCounts**: (not produced in mode *alignment-free*) This would contain the gene-level counts (output of `featureCounts <http://bioinf.wehi.edu.au/featureCounts/>`__), on the filtered GTF files, that can be used for differential expression analysis.

* **bamCoverage**: (not produced in mode *alignment-free*) This would contain the bigWigs produced by deepTools `bamCoverage <https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html>`__ . Files with suffix ``.coverage.bw`` are raw coverage files, while the files with suffix ``RPKM.bw`` are `RPKM-normalized <https://www.biostars.org/p/273537/>`__ coverage files.

* **deepTools_QC**: (produced in the mode *deepTools_QC*) This contains the quality checks specific for RNA-seq, performed via deepTools. The output folders are names after various deepTools functions and the outputs are explained under `deepTools documentation <deeptools.readthedocs.io>`__. In short, they show the insert size distribution(**bamPEFragmentSize**), mapping statistics (**estimateReadFiltering**), sample-to-sample correlations and PCA (**multiBigwigSummary, plotCorrelation, plotPCA**), and read enrichment on various genic features (**plotEnrichment**)

* **DESeq2_[sampleSheet]/DESeq2_Salmon_[sampleSheet]**: (produced in the modes *alignment* or *alignment-free*, only if a sample-sheet is provided.) This folder contains the HTML result report **DESeq2_report.html**, the annotated output file from DESeq2 (**DEseq_basic_DEresults.tsv**) and normalized counts for all samples, produced via DEseq2 (**DEseq_basic_counts_DESeq2.normalized.tsv**) as well as Rdata file(**DEseq_basic_DESeq.Rdata**) with the R objects ``dds <- DESeq2::DESeq(dds)`` and ``ddr <- DDESeq2::results(dds,alpha = fdr)``. **DESeq2_[sampleSheet]** uses gene counts from ``featureCounts/counts.tsv``, whereas **DESeq2_Salmon_[sampleSheet]** uses transcript counts from ``Salmon/counts.tsv`` that are merged via tximport in R.

* **Salmon**: (produced in mode *alignment-free*) This folder contains transcript-level (``counts.tsv``)and gene-level (``counts.genes.tsv``) counts estimated by the tool `Salmon <https://salmon.readthedocs.io/en/latest/salmon.html>`__ .

* **sleuth_Salmon_[sampleSheet]** (produced in mode *alignment-free*, only if a sample-sheet is provided) This folder contains a transcript-level differential expression output produced using the tool `Sleuth <https://pachterlab.github.io/sleuth/about>`__ .


Command line options
--------------------

.. argparse::
    :func: parse_args
    :filename: ../snakePipes/workflows/RNA-seq/RNA-seq
    :prog: RNA-seq
    :nodefault:
