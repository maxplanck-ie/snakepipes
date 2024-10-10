snakePipes News
===============

snakePipes 3.0.1
________________

 * installation test for python 3.13
 * bulk mode in makePairs wf

snakePipes 3.0.0
________________

* clusteryaml omitted for profiles
* All workflows have '-' removed from their name
* toml file installation
* makePairs mode introduced

snakePipes 2.9.0
________________

* added SEACR peaks qc
* added external bed functionality for differential binding in ChIP-seq and ATAC-seq workflows
* added "allelic-counting" mode to mRNA-seq, allowing to count reads and run DGE from allelic bam files split e.g. by whatshap
* added support for custom model formula to mRNA-seq workflow
* fixed copyfile command for sampleSheet
* removed deprecated --force argument from mamba commands
* fixes #998
* fixes #997
* fixes #996
* fixes #994
* fixes #1000
* fixes #1001

snakePipes 2.8.1
----------------
* Boosted versions on shared_env, as python 3.7 and multiqc no longer work together.
* WGBS doesn't draw a PCA on datasets with one or two samples.
* HiCExplorer version boosted in HiC mode


snakePipes 2.8.0
----------------
* Moved SalmonIndex creation from mRNAseq and scRNAseq to createIndices.
* Changed the behaviour of snakePipes createEnvs - it is no longer possible to set condaEnvDir with this function. It is required to set it with snakePipes config beforhand, instead. To ingore what's in the defaults.yaml and overwrite the condaEnvDir value with default system conda prefix, use '--autodetectCondaEnvDir'.
* Snakemake options in the defaults.yaml are now an empty string. The required arguments '--use-conda --conda-prefix' have been directly added to the command string. condaEnvDir is parsed from defaults.yaml, requiring running snakePipes config first.
* Added a 'three-prime-seq' mode to mRNAseq (David Koppstein and Katarzyna Sikora).
* Added DESeq2 run on PAS clusters to the 'three-prime-seq' mode of mRNAseq.
* Added support for multiple comparison groups to ChIPseq and ATAC-seq.
* Added SEACR as an optional peak caller to ChIPseq.
* Fixes #819
* Fixes #947
* Fixes #945
* Fixes #941
* Fixes #951
* fastq files are checked for validity
* an 'on success' file is touched in the output directory when a workflow is finished successfully
* fuzzywuzzy deprecated in favor for thefuzz


snakePipes 2.7.3
----------------
* Fixes #884 by creating an additional conda env for DSS. Build tests are now run with strict channel priority.

snakePipes 2.7.2
----------------
* STAR version has been updated to 2.7.10b. 2.7.10a was returning segmentation fault on MAC.
* STAR command has been updated. Now, STAR itself offers a command line option for processing input files.
* Put a cap on python version for the deeptools env. The current version of deeptools is not supporting the newer python versions and some tools fail.
* Update default condaDir.
* The filter_gtf function has become a bit more versatile. GTF files that include delimiters (';') in e.g. a description field are now allowed. Gene names are also allowed to have symbols now. Lastly, GTF files that have xRNA instead of transcript as a feature in column 3 can also be parsed.

snakePipes 2.7.1
----------------
* Capped snakemake version to < 7.19.0.


snakePipes 2.7.0
----------------

* Added the allelic version of Salmon-based transcript quantitation to mRNA-seq workflow. Will be run if *both* 'allelic-mapping' and 'alignment-free' modes are specified. An allelic version of sleuth will be run, if sample sheet is provided, as well as DESeq2 on allelic Salmon counts.


snakePipes 2.6.1
----------------

* Capped tabulate version as 0.9.0 breaks snakemake


snakePipes 2.6.0
----------------

* Added apeglm2 logFC shrinkage to allelic DESeq2 results.
* Added bwa-mem2 as an optional aligner to DNA-mapping and HiC, as well as as an additional tool to createIndices.
* Added bwameth2 as an optional tool to createIndices - this will create a bwameth index with bwa-mem2.
* Added bwameth2 as an optional aligner to WGBS - this will run bwameth with bwa-mem2 underneath.
* Updated software versions in environment yamls.
* Updated organism yamls.
* Updated CSAW output.
* Fixed a couple of issues in the ATAC-seq workflow after sofware versions update.
* Fixed genome size conversion to string.


snakePipes 2.5.4
----------------

* Fixed a number of minor GitHub issues: #791, #816, #807, #789, #783, #768, #827.
* Fixed misleading rule name for bamcoverage in atac-seq.
* Fixed conda env building on microsoft azure.
* Fixed CSAW report for ChIP-seq.


snakePipes 2.5.3
----------------

* Switched to mamba by default at ``createEnvs`` function.


snakePipes 2.5.2
----------------

* snakemake version >= 6.2.1
* python version >= 3.7
* added snakeMake in readme
* bug fix 777, 781
* multiqc version = 1.10.1 to be able to report QC for HiC pipeline
* added snakemake catalog yaml


snakePipes 2.5.1
----------------

* Updated Bowtie2 parameters for the cut and tag data
* Updated multibamSummary in ChIPSeq pipeline for data with spike-in
* Uncommented the BamCompare for ChIPseq pipeline with spike-in
* set a max thread for each rule
* External PR #764: always check for >0 when generating heatmap for the differential analysis done by deseq2

snakePipes 2.5.0
----------------

* Added tbb = 2020.2 to dna_mapping, creatIndices and rnae_seq yaml file
* Added cut and tag parameters to macs2 and bowtei2. The parameters have been used in Kaya-Okur et al. 2019 and can be called by using --cut_n_tag.
* Updated azure tests. python37 create envs constantly failing due its long run time. This test is now split into smaller chunks. set_macos is removed since it was completely redundant with the set.yaml


snakePipes 2.4.3
----------------

* Fixed noncoding-RNA-seq workflow without sample sheet.
* Updated links to prebuilt indices on zenodo.
* Fixed read length estimation for rMats.
* Cutadapt is set as default read trimming program for both noncoding-RNA-seq and mRNA-seq workflows.


snakePipes 2.4.2
----------------

* Deeptools coverage RPKM in mRNA-seq and noncoding-RNA-seq worflows now respects blacklist and ingoreForNorm arguments.
* In mRNA-seq and noncoding-RNA-seq workflow, deeptools qc will now also output DESeq2 size factor-normalized bigwig files.
* Fixed conda env for WGBS.
* Fixed control group ordering in split sample sheets in mRNA-seq and other workflows.
* Removed rule moving bams from allelic mRNA-seq and DNA-mapping workflows.

snakePipes 2.4.1
----------------

* Fixed sampleSheet splitting for multiple pairwise comparisons when group "All" is not listed.

snakePipes 2.4.0
----------------

* Added support for multiple pairwise comparisons for DESeq2, sleuth, and rMats in the mRNA-seq workflow, as well as for DESeq2 in the noncoding-RNA-seq workflow.
* Loompy from conda is now used in mode STARsolo in scRNAseq workflow.
* Added bamExt to mRNA-seq and noncoding-RNA-seq commandline arguments.
* Added multi-thread support to rMats in mRNA-seq workflow.
* Fixed deepTools GC bias command with SE reads.
* Bumped HiC explorer version.
* Fixed STARsoloCoords for Custom kit.


snakePipes 2.3.1
----------------

* Fixed aligner options for bwa in DNA-mapping.
* Fixed allelic mode for single end reads .
* Bumped hiC explorer version in HiC.


snakePipes 2.3.0
----------------

* Deprecated mode Gruen in scRNAseq.
* scRNAseq mode Alevin now outputs spliced/unspliced counts for RNA velocity estimation based on Soneson et al.  2020, bioRxiv https://doi.org/10.1101/2020.03.13.990069 .
* Fixed "external_gene_name" and "Status" columns in DESeq2 html report.
* Removed warning when sample names start with a number.


snakePipes 2.2.3
----------------

* Genrich will now run if sampleSheet without replicates is provided.
* Updated zenodo link to mouse genome GRCm38/mm10 .
* Fixed start coordinates in Filtered results bed from CSAW.


snakePipes 2.2.2
----------------

* Fix DAG inconsistencies for ChIP-seq and ATAC-seq ran fromBAM and from -d.
* DESeq2 Rmd file is not deleted anymore in noncoding-RNAseq.
* Fixed labels in deepTools commands.
* Allele_info is now boolean.


snakePipes 2.2.1
----------------

* Fix a bug in DAG for ChIPseq allelic with CSAW.
* Fixed deepTools qc DAG for ChIPseq with spikein.
* Added DAG test for allelic ChIPseq.
* Fixed a bug with deepTools QC for allelic mRNAseq.


snakePipes 2.2.0
----------------
* Added Alevin mode in scRNA workflow
* Added a new conda environment using to call AlevinQC.
* Added filtering of empty drops with Dropletutils to scRNA-seq mode STARsolo
* Added spikein normalization to ChIPseq workflow
* Added hybrid genome creation to createIndices
* Added STARsolo report for all samples to STARsolo output folder
* FASTQ1 and FASTQ2 are not localrules anymore due to buggy logging
* Included optional differential splicing analysis using rmats within mRNA-seq workflow
* Symlinks in the output path are relative
* Increased BBmap version
* Increased STAR version to 2.7.4a in scRNAseq, noncoding-RNA-seq and mRNA-seq workflows
* Fixed snakemake version at 5.18.0 due to a bug in DAG handling
* Minor changes to shared FastQC and multiQC rule with regards to scRNA-seq workflow.
* Fixed issue with missing input for running the DNA-mapping Snakefile
* Fixed rule TrimGalore for single end reads
* deepTools heatmaps for differentially bound regions are now ordered by sample sheet condition
* Genrich is now run on namesorted bams
* Workflow help message now points to example sampleSheet on GitHub
* organismsDir can now be updated with snakePipes config mode "recycle"

.. note::
   Please be aware that this version requires regeneration of STAR indices!

snakePipes 2.1.2
----------------
* small bug fix: SE mode in noncoding-RNA-seq pipeline

snakePipes 2.1.1
----------------
* small bug fix: a typo in atac-seq pipeline

snakePipes 2.1.0
----------------

 * Snakemake version is bumped to 5.13.0
 * Updated docs on running single snakefiles
 * Added user-input target regions and freetext parameters to differential methylation analysis with metilene
 * Added PCA to metilene report in WGBS
 * Added Genrich support for SE data
 * Edited symlinking rules to `ln -s` or python
 * TMPDIR is now passed at rule-level to the shell
 * Added logs in a couple of places
 * Added `--skipBamQC` to WGBS to be included with `--fromBAM` to suppress recalculation of QC metrics on the bam file
 * Added tempDir check to snakePipes info
 * Added `--oldConfig` and `--configMode` options to snakePipes config that allow passing a copy of an existing pre-configured config file instead of passing the single paths. Previous mode can be used with `--configMode manual` (default), the new mode with `--configMode recycle`.
 * Updated histoneHMM version to 1.8. Changed number formatting in histoneHMM output from scientific to general.
 * Small fixes in DESeq2 report for noncoding-RNA-seq, WGBS reports
 * Fixed `--verbose` in WGBS
 * Fixed an important bug in differential binding analysis with CSAW (mismatch between sampleSheet rownames and countdata colnames).


snakePipes 2.0.2
----------------

 * DAG print is now moved to _after_ workflow run execution such that any error messages from e.g. input file evaluation do not interfere with the DAG and are visible to the user.
 * Fixed fastqc for --forBAM .
 * Fixed DESeq2 report failure with just 1 DEG.
 * Updated links to test data and commands on zenodo in the docs.
 * SampleSheet check now explicitly checks for tab-delimited header.
 * Fixed metilene groups, as well methylation density plots in WGBS.

snakePipes 2.0.1
----------------

 * Fixed a bug in `snakePipes config` that caused the `toolsVersion` variable to be removed from `defaults.yaml`. This is likely related to issue #579.

snakePipes 2.0.0
----------------

 * Added a noncoding-RNA-seq workflow and renamed RNA-seq to mRNA-seq for clarity. The noncoding workflow will also quantify protein coding genes, but its primary use is analyzing repeat expression.
 * In order to use the noncoding-RNA-seq workflow organism YAML files must now include a `rmsk_file` entry.
 * Fixed STAR on CIFS mounted VFAT file systems (issue #537).
 * Added mode STARsolo to scRNAseq. This mode is now default.
 * Added log fold change shrinkage with "apeglm" to DESeq2 basic in the mRNAseq workflow. Two versions of results tables (with and without shrinkage) are now written to the DESeq2 output folder.
 * Added Genrich as peakCaller option to ChIPseq and ATACseq.
 * Added HMMRATAC as peakCaller option to ATACseq.
 * ATAC-seq short bam (filtered for short fragments) is now stored in a separate folder.

.. note::
   Please be aware that this version requires regeneration of STAR indices!

snakePipes 1.3.2
----------------

 * Fixed missing multiQC input in allelic RNAseq
 * Added sample check to those workflows that were missing it.

snakePipes 1.3.1
----------------

 * Support for snakeMake 5.7.0

snakePipes 1.3.0
----------------

 * Overhauled WGBS pipeline
 * Standardized options to be camelCase
 * Further standardized options between pipelines
 * UMI handling is now available in most pipelines
 * The ``--fromBAM`` option is now available and documented
 * Users can now change the read number indicator ("_R1" and "_R2" by default) as well as the fastq file extension on the command line.
 * Added the preprocessing pipeline, prevented python packages in users' home directories from inadvertently being used.
 * Added a ``snakePipes config`` command that can be used in lieu of editing defaults.yaml

snakePipes published
--------------------
snakePipes was published: https://www.ncbi.nlm.nih.gov/pubmed/31134269

snakePipes 1.2.3
----------------

 * Updated citation for snakePipes
 * Fixed replicate check for samples with trailing spaces in names
 * Fixed input filtering in CSAW
 * Several allele-specific RNAseq fixes
 * ATACseq peakQC is now run on fragment-size filtered bam
 * Fixed Salmon output (Number of Reads output in "prefix_counts.tsv" files and file naming)
 * Fixed CSAW QC plot error with single end reads
 * Updated histone HMM environment to a working conda version
 * Salmon_wasabi is now a localrule


snakePipes 1.2.2
----------------

 * Fixed a bug in the ATAC-seq environment where GenomeInfoDbData was missing.
 * Also an occasional issue with CSAW


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
   * User can select a metric to maximize during cell filtering (cellFilterMetric, default: gene_universe).
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
