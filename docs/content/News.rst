snakePipes News
===============

snakePipes 1.1.0
----------------

**26.11.18**

 * A wide number of bug fixes to scRNA-seq and other pipelines. In particular, many memory limits were updated.
 * An optional email can be sent upon pipeline completion.
 * The RNA-seq pipeline can now produce a fuller report upon completion if you are performing differential expression.
 * Sample merging in HiC works properly.
 * GTF files are now handled more generically, which means that they no longer need to have "_gencode" and "_ensembl" in their path.
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

snakePipes 1.0.0 (king cobra) released
--------------------------------------

**9.10.2018**

First stable version of snakePipes has been released with various feature improvements. You can download it `from GitHub<https://github.com/maxplanck-ie/snakepipes/releases/tag/1.0.0>`__

snakePipes preprint released
----------------------------

We relased the preprint of snakePipes describing the implementation and usefullness of this tool in integrative epigenomics analysis. `Read the preprint on bioRxiv <https://www.biorxiv.org/content/early/2018/09/04/407312>`__
