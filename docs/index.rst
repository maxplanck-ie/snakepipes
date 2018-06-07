snakePipes
=============

snakePipes are pipelines built using snakemake for the analysis of various sequencing datasets.

Below is the list of pipelines available in snakePipes
-------------------------------------------------------


=============================== ===========================================================================================
Pipeline                            Description
=============================== ===========================================================================================
:ref:`DNA-mapping`              Basic DNA mapping using bowtie2, filter mapped files, QC and create coverage plots
:ref:`ChIP-Seq`                 Use the DNA mapping output and run ChIP/Input normalization and peak calling
:ref:`ATAC-seq`                 Use the DNA mapping output and detect open chromatin regions from ATAC-seq
:ref:`RNA-Seq`                  RNA-Seq workflow : From mapping to differential expression using DEseq2
:ref:`scRNA-Seq`                Single-cell RNA-Seq workflow : From mapping to differential expression
:ref:`HiC`                      Hi-C analysis workflow, from mapping to TAD calling
:ref:`WGBS`                     WGBS analysis workflow, from mapping to DMR calling and differential methylation analysis
=============================== ===========================================================================================

Quick start
----------------

The desing of all workflows, as well as configuration is the same. Here's an example with ChIP-seq data.

A **typical ChIP-seq analysis of human samples** starting from paired-end FASTQ files in the directory `input-dir`:

.. code:: bash

    $ ls /path/to/input-dir/
    my_H3K27ac_sample_R1.fastq.gz  my_H3K27me3_sample_R1.fastq.gz  my_Input_sample_R1.fastq.gz
    my_H3K27ac_sample_R2.fastq.gz  my_H3K27me3_sample_R2.fastq.gz  my_Input_sample_R2.fastq.gz
    $
    $ snakemake_workflows/DNA-mapping \
          -i /path/to/input-dir -o /path/to/output-dir --dedup hs37d5 && \
      snakemake_workflows/ChIP-seq chip-seq.config.yaml \
          -d /path/to/outputdir hs37d5

Here, hs37d5 is the name of the genome. The yaml file corresponding to this genome should exist under `/shared/organisms/hs37d5.yaml`.
This yaml file should have paths to the required genome fasta, index, GTF and other annotations (see **Genome configuration file** below).

All individual jobs of the workflow will be submitted to the Grid engine using the command specified under /shared/cluster.yaml.
To run the workflow locally, use the parameter `--local` for local mode and the parameter `-j 48` to specify the maximal
number of used CPU threads (here: 48) or concurrent running Slurm jobs (actual used threads are defined in each rule).

A **configuration file is required for the ChIP-seq workflow** and should adhere to the following style :
**IMPORTANT: Use only whitespace, but NO TAB indentation in this file:**

.. code:: bash

    $ cat chip-seq.config.yaml

.. parsed-literal::

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


Genome configuration file
----------------------------

Any organism can be supported in the workflows by adding a genome configuration file `my_organism.yaml` in the
following style to the `snakemake_workflows` directory:

.. code:: bash

    $ cat snakemake_workflows/shared/organisms/hs37d5.yaml

.. parsed-literal::

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


.. note:: If no blacklist regions are available for your organism of interest, leave `blacklist_bed:` empty


Contents:
------------

.. toctree::
   :maxdepth: 2

   content/setting_up.rst
   content/workflows/DNA-mapping.rst
   content/workflows/ChIP-seq.rst
   content/workflows/RNA-seq.rst
   content/workflows/scRNA-seq.rst
   content/workflows/HiC.rst
   ChangeLog.rst

Citation
---------

.. image:: content/images/logo_mpi-ie.jpg

This tool suite is developed by the `Bioinformatics Unit <http://www.ie-freiburg.mpg.de/bioinformaticsfac>`_
at the `Max Planck Institute for Immunobiology and Epigenetics <http://www.ie-freiburg.mpg.de/>`_, Freiburg.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
