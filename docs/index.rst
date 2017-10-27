snakePipes
=============

snakePipes are pipelines built using snakemake for the analysis of various sequencing datasets.

The following is the list of pipelines available in snakePipes
---------------------------------------------------------------


=============================== ===========================================================================================
Pipeline                            Description
=============================== ===========================================================================================
:ref:`DNA-mapping`              Basic DNA mapping using bowtie2, filter mapped files, QC and create coverage plots
:ref:`ChIP-Seq`                 Use the DNA mapping output and run ChIP/Input normalization and peak calling
:ref:`RNA-Seq`                  RNA-Seq workflow : From mapping to differential expression using DEseq2
:ref:`scRNA-Seq`                Single-cell RNA-Seq workflow : From mapping to differential expression
=============================== ===========================================================================================


Contents:
------------

.. toctree::
   :maxdepth: 2

   content/workflows/DNA-mapping.rst
   content/workflows/ChIP-seq.rst
   content/workflows/RNA-seq.rst
   content/workflows/scRNA-seq.rst

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
