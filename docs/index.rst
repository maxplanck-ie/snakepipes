snakePipes
==========

.. image:: content/images/snakePipes_small.png
   :scale: 20 %
   :height: 200px
   :width: 200 px

snakePipes are pipelines built using `snakemake <snakemake.readthedocs.io>`__ and *python* for the analysis of epigenomic datasets.


Below is the list of pipelines available in snakePipes
------------------------------------------------------

=============================== ===============================================================================================================
Pipeline                            Description
=============================== ===============================================================================================================
:ref:`createIndices`            Create indices for an organism for further use within snakePipes
:ref:`DNA-mapping`              Basic DNA mapping using bowtie2, filter mapped files, QC and create coverage plots
:ref:`ChIP-Seq`                 Use the DNA mapping output and run ChIP/Input normalization and peak calling
:ref:`ATAC-seq`                 Use the DNA mapping output and detect open chromatin regions for ATAC-seq data
:ref:`HiC`                      Hi-C analysis workflow, from mapping to TAD calling
:ref:`noncoding-RNA-Seq`        noncoding-RNA-Seq workflow : From mapping to differential expression of genes and repeat elements using DESeq2
:ref:`mRNA-Seq`                 RNA-Seq workflow : From mapping to differential expression using DESeq2
:ref:`scRNA-Seq`                Single-cell RNA-Seq (CEL-Seq2) workflow : From mapping to differential expression
:ref:`WGBS`                     Whole-genome Bisulfite-Seq analysis workflow, from mapping to DMR calling and differential methylation analysis
:ref:`preprocessing`            Merging technical replicates (e.g., across lanes), removing optical duplicates, running FastQC
=============================== ===============================================================================================================

Quick start 
-----------

* Assuming you have *python3* with *conda*, install snakePipes with:

.. code:: bash

    conda create -n snakePipes -c mpi-ie -c bioconda -c conda-forge snakePipes==2.1.1

* You can update snakePipes to the latest version available on conda with:

.. code:: bash

    conda update -n snakePipes -c mpi-ie -c bioconda -c conda-forge --prune snakePipes

snakePipes is going to move to mamba in the long-term perpective.

* Download genome fasta and annotations for an your organism, and build indexes, Check in :ref:`createIndices`

* Configure snakePipes with paths to organism and cluster configs on your system using snakePipes config. For detailed information, run:

.. code:: bash

    snakePipes config --help

.. note:: If you have a copy of a `shared/defaults.yaml` with the necessary paths configured (i.e. from a previous installation), you can pass it to snakePipes config with `--oldConfig` and `--configMode recycle` instead of providing all the paths manually again. Config keys have to match for this to work.

* Download example fastq files for the human genome `here <https://zenodo.org/record/3707259>`_

* Execute the DNA-mapping pipeline using the example **command.sh** in the test data directory.


Running your own analysis
-------------------------

For a detail introduction to setting up snakePipes from scratch, please visit :doc:`content/setting_up`

For each organism of interest, snakePipes requires fasta files, genome indexes and annotation files.
Paths to these files are specified in the organism/<name>.yaml files. After installation, the location
of these files could be revealed by the following command:

.. code:: bash

    snakePipes info

You could either modify the existing files (add your own paths), or add a new file there. See more detail in
:doc:`content/running_snakePipes`

snakePipes could either be executed locally, or on any snakemake-supported cluster infrastructure. See details
for setting up the cluster command in :doc:`content/running_snakePipes`

Citation
-------------

If you adopt/run snakePipes for your analysis, cite it as follows :

Bhardwaj V, Heyne S, Sikora K, Rabbani L, Rauer M, Kilpert F, Richter AS, Ryan DP, Manke T. **snakePipes: facilitating flexible, scalable and integrative epigenomic analysis.** Bioinformatics. 2019 May 27. pii: btz436. doi:

10.1093/bioinformatics/btz436. [Epub ahead of print] PubMed PMID: 31134269. https://www.ncbi.nlm.nih.gov/pubmed/31134269


.. image:: content/images/logo_mpi-ie.jpg

This tool suite is developed by the `Bioinformatics Unit <http://www.ie-freiburg.mpg.de/bioinformaticsfac>`_
at the `Max Planck Institute for Immunobiology and Epigenetics <http://www.ie-freiburg.mpg.de/>`_, Freiburg.


Help and Support
----------------

For query/questions regarding snakePipes, please write on biostars with the tag **#snakePipes**

For feature requests or bug reports, please open an issue on `our GitHub Repository <https://github.com/maxplanck-ie/snakepipes>`__.


Contents:
---------

.. toctree::
   :maxdepth: 2

   content/setting_up.rst
   content/running_snakePipes.rst
   content/advanced_usage.rst
   content/workflows/createIndices.rst
   content/workflows/DNA-mapping.rst
   content/workflows/ChIP-seq.rst
   content/workflows/ATAC-seq.rst
   content/workflows/HiC.rst
   content/workflows/preprocessing.rst
   content/workflows/mRNA-seq.rst
   content/workflows/noncoding-RNA-seq.rst
   content/workflows/scRNA-seq.rst
   content/workflows/WGBS.rst
   content/News.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
