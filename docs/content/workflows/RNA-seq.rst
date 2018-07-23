.. _RNA-seq:

RNA-seq
=======

What it does
------------

The snakePipes RNA-seq workflow allows users to process their single or paired-end
RNA-Seq fastq files upto the point of gene/transcript-counts and differential expression.
It also allows full allele-specific RNA-Seq analysis (up to allele-specific
differential expression) using the *allelic-mapping* mode.

.. image:: ../images/RNAseq_pipeline.png

Input requirements and output
-----------------------------

The only requirement is a directory of gzipped fastq files. Additionally for
differential expression analysis, a tab-separated sample information file could be provided.


Analysis modes
--------------

Following analysis (*modes*) are possible using the RNA-seq workflow:

"alignment"
~~~~~~~~~~~

In this mode,
the pipeline uses one of the selected aligners to create BAM files, followed by
gene-level quantification using **featurecounts**. Gene-level differential expression
analysis is then performed using **DESeq2**.

"allelic-mapping"
~~~~~~~~~~~~~~~~~

**allelic-mapping** mode follows a similar process as the "mapping" mode, however the
alignment performed on an allele-masked genome, followed by allele-specific splitting
of mapped files. Gene-level quantification is performed for each allele using *featureCounts*.
Allele-specific, gene-level differential expression analysis is then performed using **DESeq2**.

.. note:: **allelic-mapping** mode is mutually exclusive with **mapping** mode

"alignment-free"
~~~~~~~~~~~~~~~~

In this mode,
the pipeline uses **salmon** to perform transcript-level expression quantification. This mode
performs both transcript-level differential expression (using **Sleuth**), and gene-level
differential expression (using **wasabi**, followed by **DESeq2**).

"deepTools_qc"
~~~~~~~~~~~~~~

The pipeline provides multiple quality controls through deepTools, which can be triggered
using the *deepTools_qc* mode. It's a very useful add-on with any of the other modes.


Command line options
--------------------

.. argparse::
    :func: parse_args
    :filename: ../snakePipes/workflows/RNA-seq/RNA-seq
    :prog: RNA-seq
    :nodefault:
