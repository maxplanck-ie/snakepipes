.. _HiC:

HiC
============

Input requirements and outputs:
-------------------------------------------
This pipeline requires paired-end reads fastq files as input to build a contact matrix and to call TADs.


.. argparse::
   :filename: ../workflows/HiC/HiC
   :func: parse_args
   :prog: HiC


Structure of output directory
-------------------------------

The HiC pipeline will generate output of the following structure:

.. code:: bash

    $ tree -d -L 2 output-dir/


.. parsed-literal::

 output-dir/
    |-- BWA
    |-- FASTQ
    |-- HiC_matrices
    |     |-- logs
    |     |-- QCplots
    |-- HiC_matrices_corrected
    |     |-- logs
    |-- TADs
          |--logs


