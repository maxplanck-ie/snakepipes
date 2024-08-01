.. _scRNAseq:

scRNAseq
=========

What it does
------------

The scRNAseq pipeline is intended to process UMI-based data, expecting the cell barcode and umi in Read1, and the cDNA sequence in Read2. The workflow has predefined settings for CelSeq2 and 10x data, but can be extended to custom protocols.

There are currently two analysis modes available:
- "STARsolo" which uses STAR solo for mapping and quantitation.
- "Alevin" based on Salmon for generating the count matrix.

.. note:: Mode "Gruen" has been deprecated.

The general procedure for mode "STARsolo" involves:

1. moving cell barcodes and UMIs from read 1 into the CB and UMI tags of read 2 during mapping (STARsolo),
2. quantification of genic read counts at the single cell level (STARsolo),
3. quantification of reads supporting spliced and unspliced transcripts in each cell (velocyto) - unless this has been disabled with --skipVelocyto
4. generation of seurat objects for genic counts.

UMIs in the read headers are used to avoid counting PCR duplicates. A number of bigWig and QC plots (e.g., from ``plotEnrichment``) are generated as well.

Mode "Alevin" involves:

1. Generation of a salmon index used for mapping.
2. Mapping and generation of a readcount matrix.
3. Estimation of uncertainty of gene counts using bootstrap method implemented in Salmon Alevin.
4. General QC of the Alevin run using the AlevinQC R package.
5. Quantification of "spliced" and "unspliced" read counts in each cell with Alevin - unless this has been disabled with --skipVelocyto . This analysis is derived from the code underlying Soneson et al.  2020, bioRxiv https://doi.org/10.1101/2020.03.13.990069. 

.. image:: ../images/scRNAseq_pipeline.png


Mode STARsolo
-------------

With current settings, this mode should work with any UMI-based protocol that stores UMI and CB in read 1, each in one chunk. 
The mode comes with four presets that can be passed to the ' --myKit ' argument: CellSeq192, CellSeq384, 10xV2, 10xV3. Choosing a preset will select a corresponding barcode whitelist file as well as cell barcode and umi length and positions to be used. Choosing the Custom preset allows the user to run the workflow providing own barcode whitelist and CB/UMI positions and lengths. CellSeq384 is the current default preset.

In this mode, STARsolo is used to map, UMI-deduplicate and count reads. Importantly, read 1 is expected to carry the UMI and the cell barcode, while read 2 is expected to carry the cDNA sequence. Default positions of UMI and CB in read 1 are specified, as well as their respective lengths. If your setup is different from the available presets, change it via the --STARsoloCoords commandline argument or in the defaults.yaml dictionary, in addition to providing --myKit Custom argument.

In the STARsolo folder, bam files are stored, along with 10x-format count matrices and log files summarizing barcode detection and UMI-deduplication.
Bam files have the UB and CB tags set.

Deeptools QC is run on these bam files.

Before running velocyto, bam files from STARsolo are filtered to remove unmapped reads as well as reads with an empty CB tag and then cell-sorted by the CB tag.
In the VelocytoCounts folder, loom files with counts of spliced, unspliced and ambiguous reads are stored. A merged loom file containing counts for all samples together can be found in the VelocytoCounts_merged folder. As Velocyto tends to consume a lot of memory and result in long runtimes with cell numbers in ~10^5, it can be disabled with --skipVelocyto.


Input requirements
------------------

The primary input requirement is a directory of paired-end fastq files. For the STAR solo mode, a barcode whitelist is required, as well as specification of UMI and CB positions and length, if different from default or available presets.

Barcode whitelist
~~~~~~~~~~~~~~~~~

Required for the STARsolo mode. The expected format is a one-column txt file with barcodes the user wishes to retain. Default is a whitelist file for CellSeq2 384 barcodes, provided with the pipeline. If 'myKit' is changed to another available preset, the corresponding barcode whitelist provided with the pipeline will be used.


Configuration file
~~~~~~~~~~~~~~~~~~

The default configuration file is listed below and can be found in ``snakePipes/workflows/scRNAseq/defaults.yaml``::

    pipeline: scrnaseq
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
    ##Analysis mode
    mode: STARsolo
    ## Number of reads to downsample from each FASTQ file
    downsample:
    ## Options for trimming
    trim: False
    trimmer: cutadapt
    trimmerOptions: -a A{'30'}
    ## --twopassMode Basic is not compatible with --outStd in all STAR versions
    alignerOptions: ""
    ## further options
    filterGTF: "-v -P 'decay|pseudogene' "
    cellBarcodeFile:
    cellBarcodePattern: "NNNNNNXXXXXX"
    splitLib: False
    cellNames:
    ##mode STARsolo options
    myKit: CellSeq384
    BCwhiteList:
    STARsoloCoords: ["1","7","8","7"]
    skipVelocyto: False
    ##mode Alevin options
    alevinLibraryType: "ISR"
    prepProtocol: "chromiumV3"
    expectCells: 
    readLengthFrx: 0.2
    #generic options
    libraryType: 1
    bwBinSize: 10
    verbose: False
    plotFormat: pdf
    dnaContam: False
    ## Parameters for th statistical analysis
    cellFilterMetric: gene_universe
    #Option to skip RaceID to save time
    skipRaceID: False
    #umi_tools options:
    UMIBarcode: False
    bcPattern: NNNNCCCCCCCCC #default: 4 base umi barcode, 9 base cell barcode (eg. RELACS barcode)
    UMIDedup: False
    UMIDedupSep: "_"
    UMIDedupOpts: --pairedUMIDedupOpts: --paired


Pseudogene filter
~~~~~~~~~~~~~~~~~

As default, transcripts or genes that contain that are related to biotypes like 'pseudogene' or 'decay' are filtered out before tag counting (see
:code:`--filterGTF` default).
Here we assume you provide eg. a gencode or ensemble annotation file (via genes_gtf in the organism configuration yaml) that contains this information.

Library Type
~~~~~~~~~~~~

The CEL-seq2 protocol produces reads where read 2 maps in sense direction (:code:`libraryType: 1`).


Fraction of read length required to overlap the intron
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In mode Alevin, the fraction of read length required to overlap the intron in order to be counted as "unspliced" is set to 0.2 (i.e. 20%) by default. This corresponds to 10nt in a 50nt-long read, or to 20nt in a 100nt-long read. The user is encouraged to modify this value as deemed appropriate via the ``--readLengthFrx`` commandline argument.
In practice, this variable affects the length of the exon sequence flank added to the intron sequence to generate reference sequences for Salmon Alevin. Exon sequence flank length is set to one minus 'readLengthFrx' of read length.



Output structure
----------------

The following will be produced in the output directory when the workflow is run in mode STARsolo::

    analysis/
    ├── scRNAseq_run-1.log
    ├── multiQC
    ├── deepTools_qc
    ├── cluster_logs
    ├── bamCoverage
    ├── Sambamba
    ├── filtered_bam
    ├── STARsolo
    ├── Seurat
    ├── Annotation
    ├── FastQC
    ├── originalFASTQ
    ├── scRNAseq_tools.txt
    ├── scRNAseq.cluster_config.yaml
    ├── scRNAseq.config.yaml
    └── scRNAseq_organism.yaml

 - The **VelocytoCounts** directory contains loom files in sample subdirectories.
 - The **VelocytoCounts_merged** directory containes one loom file with all samples merged.
 - The **STARsolo** directory contains bam files and 10X-format cell count matrices produced by STARsolo.
 - The **Annotation** directory contains a filtered version of your original GTF file, with pseudogenes removed by default.
 - The **bamCoverage** directory contains a bigwig track for each sample (not per cell!). This can be used eg. in IGV to check where your reads map in general.
 - The **deeptools_qc** directory contains additional QC reports and plots. The ``FASTQC`` directory can be used to verify eg. the barcode layout of read 1.


The following output structure will be produced when running in Alevin mode::

    ├── Alevin
    ├── AlevinForVelocity
    ├── Annotation
    ├── cluster_logs
    ├── FastQC
    ├── multiQC
    ├── originalFASTQ
    ├── Salmon
    ├── scRNAseq.cluster_config.yaml
    ├── scRNAseq.config.yaml
    ├── scRNAseq_organism.yaml
    ├── scRNAseq_pipeline.pdf
    ├── scRNAseq_run-1.log
    ├── scRNAseq_tools.txt
    └── SingleCellExperiment

 - The **Salmon** directory contains the generated genome index.
 - The **Alevin** directory contains the matrix files (both bootstrapped and raw) per sample in subdirectories.
 - The **multiQC** directory contains an additional alevinQC html file generated per sample.
 - The **AlevinForVelocity** directory contains the matrix files with "spliced" and "unspliced" reads per cell in subdirectories.
 - The **SingleCellExperiment** directory contains the RDS files with "SingleCellExperiment" class R objects, storing spliced/unspliced counts per cell in corresponding assays.

Understanding the outputs: mode STARsolo
----------------------------------------

- **Main result:** output folders with 10x-format count matrices can be found in sample subfolders under ``STARsolo``. The ouput consists of three files: barcodes.tsv, features.tsv, matrix.mtx. Their gzipped versions are stored in the same folder. Seurat objects from merged samples are available in the ``Seurat`` folder.

- Corresponding annotation files are: ``Annotation/genes.filtered.bed`` and ``Annotation/genes.filtered.gtf``, respectively.

- The folders ``QC_report``, ``FASTQC``, ``deeptools_qc`` and ``multiQC`` contain various QC tables and plots.

- *STARsolo* directory contain the output from genomic alignments.


Understanding the outputs: mode Alevin
--------------------------------------

- **Main result:** output folders containing the raw and boostrapped count matrices are found under the sample subfolders under ``Alevin``. The sample specific Alevin folders contain the matrices, as well as column data (barcodes) and row data (genes). Alevin spliced/unspliced counts for RNA velocity are stored as alevin matrices in the sample subfolders under ``AlevinForVelocity`` and as "SingleCellExperiment" class R objects under ``SingleCellExperiment``.

- Corresponding annotation files are: ``Annotation/genes.filtered.bed`` and ``Annotation/genes.filtered.gtf``, respectively.

- The QC plots (both from multiQC and AlevinQC) are available in the ``multiQC`` folder.


Command line options
--------------------

.. argparse::
    :func: parse_args
    :filename: ../snakePipes/workflows/scRNAseq/scRNAseq
    :prog: scRNAseq
    :nodefault:
