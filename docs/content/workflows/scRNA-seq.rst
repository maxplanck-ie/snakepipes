.. _scRNA-seq:

scRNA-seq
=========

What it does
------------

The scRNA-seq pipeline is intended to process CEL-Seq2 data, though it may be able to process some similar Drop-seq protocols. The general procedure involves

1. moving cell barcodes and UMIs from read 1 into the read headers of read 2,
2. mapping read 2,
3. quantification at the single cell level.

UMIs in the read headers are used to avoid counting PCR duplicates. A number of bigWig and QC plots (e.g., from ``plotEnrichment``) are generated as well.

.. image:: ../images/scRNAseq_pipeline.png

Input requirements
------------------

The primary input requirement is a directory of paired-end fastq files. In addition, if you do not wish to use the default list of cell-barcodes you must then supply your own.

Cell barcodes
~~~~~~~~~~~~~

The format of the cell barcodes file is shown below. Note that the default file is included in the snakePipes source code under ``snakePipes/workflows/scRNAseq``. This file is automatically used if you leave :code:`barcode_file` empty.

::

    1       AGTGTC
    2       ACCATG
    3       GAGTGA
    4       CACTCA
    5       CATGTC
    6       ACAGGA
    7       GTACCA
    8       ACAGAC
    9       ACGTTG

The default cell barcodes are 192 hexamers listed in a file with the first column a cell number and the second the barcode sequence.

Predefined cell barcodes are required right now. However it is planned to make this more generic in future workflow versions.

Configuration file
~~~~~~~~~~~~~~~~~~

The default configuration file is listed below and can be found in ``snakePipes/workflows/scRNAseq/defaults.yaml``::

    pipeline: scrna-seq
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
    ## Number of reads to downsample from each FASTQ file
    downsample:
    ## Options for trimming
    trim: False
    trimmer: cutadapt
    trimmerOptions: -a A{'30'}
    ## further options
    filter_annotation: "-v -P 'decay|pseudogene' "
    barcode_file:
    barcode_pattern: "NNNNNNXXXXXX"
    split_lib: False
    cell_names:
    libraryType: 1
    bwBinSize: 10
    verbose: False
    plotFormat: pdf
    dnaContam: False
    ## Parameters for th statistical analysis
    cell_filter_metric: gene_universe
    #Option to skip RaceID to save time
    skipRaceID: False


While some of these can be changed on the command line, you may find it useful to change ``barcode_pattern`` and ``barcode_file`` if you find that you need to change them frequently.

Barcode pattern
~~~~~~~~~~~~~~~

The scRNA-seq pipeline requires barcodes at 5' end of read 1. The default barcode_pattern takes the first 6 bases as UMI (NNNNNN) and the following 6 bases as cell barcode (XXXXXX).
If your read/barcode layout requires additional **'Don't care'** positions eg. before stretches of N one can indicate these with ``.``

Barcode file
~~~~~~~~~~~~~~~

Only specify a file if you use other than the default CEL-seq2 barcodes.


Trimming
~~~~~~~~

It is recommended to use the :code:`--trim` option as this uses cutadapt to trim remaining adapters *and* poly-A tails from read 2 (see defaults for ``--trimmerOptions``).

Pseudogene filter
~~~~~~~~~~~~~~~~~

As default, transcripts or genes that contain that are related to biotypes like 'pseudogene' or 'decay' are filtered out before tag counting (see
:code:`--filter_annotation` default).
Here we assume you provide eg. a gencode or ensemble annotation file (via genes_gtf in the organism configuration yaml) that contains this information.

Library Type
~~~~~~~~~~~~

The CEL-seq2 protocol produces reads where read 2 maps in sense direction (:code:`libraryType: 1`). After barcodes are transferred to read 2, the workflow continues in single-end mode.

Split lib
~~~~~~~~~

This option you need in case a library contains only 96 instead of 192 cells.



Output structure
----------------

The following will be produced in the output directory::

    |-- cluster_logs
    |-- Filtered_cells_RaceID
    |   `-- logs
    |-- Filtered_cells_monocle
    |    `-- logs
    |-- cellQC_test
    |-- mtab_test
    |-- QC_report
    |   `-- data
    |-- Results
    |-- Counts
    |   `-- logs
    |-- multiQC
    |   `-- multiqc_data
    |-- bamCoverage
    |   `-- logs
    |-- deepTools_qc
    |   |-- logs
    |   |-- bamPEFragmentSize
    |   |-- plotEnrichment
    |   `-- estimateReadFiltering
    |-- Sambamba
    |-- STAR_genomic
    |   |-- logs
    |   `-- GSM2668205
    |-- FastQC
    |   `-- logs
    |-- Annotation
    |-- FASTQ_barcoded
    `-- FASTQ

 - The **Annotation** directory contains a filtered version of your original GTF file, with pseudogenes removed by default.
 - The **bamCoverage** directory contains a bigwig track for each sample (not per cell!). This can be used eg. in IGV to check where your reads map in general.
 - The **Counts** directory contains 4 sets of counts: UMIs/feature/cell (.umis.txt), reads/feature/cell (.reads.txt), corrected number of UMIs/feature/cell (corrected.txt) and raw counts per cell per UMI per feature (raw_counts.txt). Of these, the values in corrected.txt should be used for further analysis and the others for quality control.
 - The **deeptools_qc** directory contains additional QC reports and plots. The ``FASTQC`` directory can be used to verify eg. the barcode layout of read 1.
 - The **QC_report** directory contains additional QC stats as tables and plots.

Understanding the outputs
-------------------------

- **Main result:** the genes per cell count table with poisson-corrected counts can be found under ``Results/all_samples.gencode_genomic.corrected_merged.csv``

- Corresponding annotation files are: ``Annotation/genes.filtered.bed`` and ``Annotation/genes.filtered.gtf``, respectively.

- The folders ``QC_report``, ``FASTQC``, ``deeptools_qc`` and ``multiQC`` contain various QC tables and plots.

- **Sambamba** and **STAR_genomic** directories contain the output file from duplicate marking and genomic alignments, respectively.

Filtered_cells_monocle
~~~~~~~~~~~~~~~~~~~~~~

The poisson-rescaled count matrix is read and converted into a monocle dataset. A range of transcript counts per cell thresholds (from 1000 to 5000 by 500) are applied to filter cells and the resulting R objects are written to minT*.mono.set.RData. For every cell filtering threshold, several metrics are collected and written to metrics.tab.txt: number of retained cells, median number of expressed genes per cell (GPC), size of the total gene universe. Plots of median GPC as well as gene universe size as functions of the cell filtering threshold are written to medGPCvsminT.downscaled.png and gene_universevsminT.downscaled.png, respectively.

The optimal cell filtering threshold for the subsequent analyses is selected as the value that results in maximizing a gene expression metric choosable from "gene_universe" (default) and "medGPC". Using gene universe tends to maximize the overall cell diversity while using median genes per cell (medGPC) maximizes the information content per cell.
Gene expression dispersions are calculated for the corresponding monocle object and the trend plot is written to mono.set.*.disp.estim.png. A first iteration of cell clustering with default settings resutls in a rho-delta plot written to mono.set.*.rho_delta.png and a tSNE plot with cell cluster colouring written to mono.set.*.tsne.auto.Cluster.png. Rho and delta are now re-evaluated and set to the 80th and the 95th percentiles of the original distributions, respectively. Cells are reclustered and the corresponding tSNE plot is written to mono.set.*.tsne.thd.Cluster.png. The monocle object containing the updated clustering information is written to minT*.mono.set.RData. It is also converted to a seurat object and the clustering information is transferred. The seurat object is saved as minT*.seuset.RData. The tSNE plot with clustering information produced with seurat is written to minT*.seuset.tSNE.png.
Top10 as well as top2 markers are calculated for each cell cluster and written to minT\*.Top10markers.txt and minT\*.Top2markers.txt, respectively. The corresponding heatmaps are written to minT\*.Top10markers.heatmap.png and minT\*.Top2markers.heatmap.png, respectively. For the top2 marker list, violin as well as feature plots are produced and saved under Top2.clu\*.violin.png and Top2.clu\*.featurePlot.png, respectively. The R session info is written to sessionInfo.txt.
Statistical procedures and results are summarized in Stats_report.html.

Filtered_cells_RaceID
~~~~~~~~~~~~~~~~~~~~~

Cell filtering, metrics collection and threshold selection are done as above only using RaceID package functions, where applicable.

Clustering is done with RaceID default settings. The fully processed RaceID object is written to sc.minT\*.RData, the tsne plot with the clustering information to sc.minT\*.tsne.clu.png.
Top 10 and top 2 markers are calculated, and the resulting plots and tables written out as above. Violin and feature plots are generated for the top2 marker list and saved to files as in the description above. Session info is written to sessionInfo.txt. Statistical procedures and results are summarized in Stats_report.html.


Example images
~~~~~~~~~~~~~~

There are a number of QC images produced by the pipeline:

.. image:: ../images/scRNAseq_UMI_plot.png

This figure plots the number of UMIs on transcripts per cell vs the number of reads aligning to transcripts. These should form a largely straight line, with the slope indicating the level of PCR duplication.

.. image:: ../images/scRNAseq_plate_abs_transcript.png

This figure shows the distribution of the number of UMIs across the single cells. Each block is a single cell and the color indicates the number of UMIs assigned to it. This is useful for flagging outlier cells.
Note: the layout corresponds to half of a 384-well plate as this is used usually for CEL-seq2. The plot can also help to see biases corresponding to the well-plate.

Command line options
--------------------

.. argparse::
    :func: parse_args
    :filename: ../snakePipes/workflows/scRNAseq/scRNAseq
    :prog: scRNAseq
    :nodefault:
