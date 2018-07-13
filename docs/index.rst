snakePipes
=============

snakePipes are pipelines built using [snakemake](snakemake.readthedocs.io) and *python* for the analysis of epigenomic datasets.

Below is the list of pipelines available in snakePipes
-------------------------------------------------------

=============================== ========================================================================================================
Pipeline                            Description
=============================== ========================================================================================================
:ref:`DNA-mapping`              Basic DNA mapping using bowtie2, filter mapped files, QC and create coverage plots
:ref:`ChIP-Seq`                 Use the DNA mapping output and run ChIP/Input normalization and peak calling
:ref:`ATAC-seq`                 Use the DNA mapping output and detect open chromatin regions for ATAC-seq data
:ref:`HiC`                      Hi-C analysis workflow, from mapping to TAD calling
:ref:`RNA-Seq`                  RNA-Seq workflow : From mapping to differential expression using DEseq2
:ref:`scRNA-Seq`                Single-cell RNA-Seq (CEL-Seq2) workflow : From mapping to differential expression
:ref:`WGBS`                     Whole-genome Bisulfite-Seq analysis workflow, from mapping to DMR calling and differential methylation analysis
=============================== ========================================================================================================

Quick start
----------------

Assuming you have *python3* with *conda*, install snakePipes with:

.. code:: bash

    conda create -n snakePipes -c mpi-ie -c bioconda -c conda-forge snakePipes

Download genome fasta, indexes and annotations for an example organism (mouse - mm10):


Download example fastq files for the mouse genome:


Execute the example RNA-Seq mapping and differential expression analysis pipeline:



Running your own analysis
-----------------------------

For a detail introduction to setting up snakePipes from scratch, please visit :doc:`/content/setting_up.rst`

For each organism of interest, snakePipes requires fasta files, genome indexes and annotation files.
Paths to these files are specified in the organism/<name>.yaml files. After installation, the location
of these files could be revealed by the following command:

.. code:: bash

    snakePipes info

You could either modify the existing files (add your own paths), or add a new file there. See more detail in
:doc:`/content/running_snakePipes.rst`

snakePipes could either be executed locally, or on any snakemake-supported cluster infrastructure. See details
for setting up the cluster command in :doc:`/content/running_snakePipes.rst`

Citation
---------

If you adopt/run snakePipes for your analysis, cite it as follows :

    snakePipes: flexible and scalable NGS analysis pipelines built using snakemake.
    The MPI-IE Bioinformatics Facility. http://doi.org/10.5281/zenodo.1146540


.. image:: content/images/logo_mpi-ie.jpg

This tool suite is developed by the `Bioinformatics Unit <http://www.ie-freiburg.mpg.de/bioinformaticsfac>`_
at the `Max Planck Institute for Immunobiology and Epigenetics <http://www.ie-freiburg.mpg.de/>`_, Freiburg.


Help and Support
-----------------

For query/questions regarding snakePipes, please write to the deepTools user group : `deeptools@googlegroups.com`

For feature requests or bug reports, please open an issue on [our GitHub Repository](https://github.com/maxplanck-ie/snakepipes)


Contents:
------------

.. toctree::
   :maxdepth: 2

   content/setting_up.rst
   content/running_snakePipes.rst
   content/workflows/DNA-mapping.rst
   content/workflows/ChIP-seq.rst
   content/workflows/ATAC-seq.rst
   content/workflows/HiC.rst
   content/workflows/RNA-seq.rst
   content/workflows/scRNA-seq.rst
   content/workflows/WGBS.rst
   ChangeLog.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
