.. _HiC:

HiC
============
The pipeline requires paired-end reads fastq files as input and generates a contact matrix. If noTAD is not mentioned it then calls the TADs.
.. argparse::
   :filename: ../workflows/HiC/HiC
   :func: parse_args
   :prog: HiC
