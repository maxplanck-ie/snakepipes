.. _HiC:

HiC
============

Input requirements and outputs:
-------------------------------------------
This pipeline requires paired-end reads fastq files as input in order to build a contact matrix and to call TADs. 
Prior to building a matrix, the pipeline maps reads against a user provided reference genome. The output of mapping step is then used for building the contact matrix.

.. argparse::
   :filename: ../workflows/HiC/HiC
   :func: parse_args
   :prog: HiC


Structure of output directory
-------------------------------

The HiC pipeline will generate output of the following structure:

.. code:: bash

    $ tree -d -L 2 outputdir
  
output directory::

    outputdir
    |--BWA
    |--FASTQ
    |--HiC_matrices
    |   |--logs
    |   |--QCplots
    |--HiC_matrices_corrected
    |   |--logs
    |--TADs
        |--logs

