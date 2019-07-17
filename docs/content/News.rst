snakePipes News
===============

snakePipes published
--------------------
snakePipes was published: https://www.ncbi.nlm.nih.gov/pubmed/31134269

snakePipes 1.2.1
----------------

 * Fixed a typo in ``createIndices``.
 * Implemented complex experimental design in RNAseq (differential gene expression), ChIP/ATACseq (differential binding).
 * Fixed an issue with ggplot2 and log transformation in RNAseq report Rmd.
 * fastqc folder is created and its content will be added to multiqc only if fastqc flag is called.
 * fastqc-trimmed folder is created and its content will be added to multiqc only if both fastqc and trim flags are called. 

snakePipes 1.2.0
----------------

 * A number of minor bug fixes across all of the pipelines
 * Updates of all tool versions and switching to R 3.5.1
 * A ``snakePipes flushOrganisms`` option was added to remove the default organism YAML files.
 * Renamed ``--notemp`` to ``--keepTemp``, which should be less confusing

snakePipes 1.1.2
----------------

 * A number of minor bug fixes and enhancements in the HiC and WGBS pipelines
 * The RNA-seq pipeline now uses samtools for sorting. This should avoid issues with STAR running out of memory during the output sorting step.
 * Increased the memory allocation for MACS2 to 8GB and bamPEFragmentSize to 3G
 * Fixed the scRNA-seq pipeline, which seems to have been broken in 1.1.1

snakePipes 1.1.1
----------------

 * Fixed some conda environments so they could all be solved in a reasonable amount of time.
 * Updated some WGBS memory limits

snakePipes 1.1.0
----------------

 * A wide number of bug fixes to scRNA-seq and other pipelines. In particular, many memory limits were updated.
 * An optional email can be sent upon pipeline completion.
 * The RNA-seq pipeline can now produce a fuller report upon completion if you are performing differential expression.
 * Sample merging in HiC works properly.
 * GTF files are now handled more generically, which means that they no longer need to have \_gencode and \_ensembl in their path.
 * WGBS:

   * Merging data from WGBS replicates is now an independent step so that dependent rules don't have to wait for successful completion of single CpG stats but can go ahead instead.
   * Filtering of differential methylation test results is now subject to two user-modifiable parameters minAbsDiff (default 0.2) and FDR (0.02) stored in defaults.yaml.
   * Metilene commandline parameters are now available in defaults.yaml. Defaults are used apart from requesting output intervals with any methylation difference (minMethDiff 0).
   * Additional diagnostic plots are generated - p value distribution before and after BH adjustment as well as a volcano plot.
   * Automatic reports are generated in every folder containing results of statistical analysis (single CpG stats, metilene DMR stats, user interval aggregate stats), as long as sample sheet is provided.
   * R sessionInfo() is now printed at the end of the statistical analysis.

 * scRNAseq:

   * An extention to the pipeline now takes the processed csv file from Results folder as input and runs cell filtering with a range of total transcript thresholds using monocle and subsequently runs clustering, produces tsne visualizations, calculates top 2 and top10 markers per cluster and produces heatmap visualizations for these using monocle/seurat. If the skipRaceID flag is set to False (default), all of the above are also executed using RaceID.
   * Stats reports were implemented for RaceID and Monocle/Seurat so that folders Filtered_cells_RaceID and Filtered_cells_monocle now contain a Stats_report.html.
   * User can select a metric to maximize during cell filtering (cell_filter_metric, default: gene_universe).
   * For calculating median GPC, RaceID counts are multiplied by the TPC threshold applied (similar to 'downscaling' in RaceID2).

 * all sample sheets now need to have a "name" and a "condition" column, that was not consistent before
 * consistent --sampleSheet [FILE] options to invoke differential analysis mode (RNA-seq, ChIP-seq, ATAC-seq), --DE/--DB were dropped

snakePipes 1.0.0 (king cobra) released
--------------------------------------

**9.10.2018**

First stable version of snakePipes has been released with various feature improvements. You can download it `from GitHub <https://github.com/maxplanck-ie/snakepipes/releases/tag/1.0.0>`__

snakePipes preprint released
----------------------------

We relased the preprint of snakePipes describing the implementation and usefullness of this tool in integrative epigenomics analysis. `Read the preprint on bioRxiv <https://www.biorxiv.org/content/early/2018/09/04/407312>`__
