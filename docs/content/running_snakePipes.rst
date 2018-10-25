Running snakePipes
==================

Pipelines under snakePipes are designed in a way such that all workflows are configured and ran in a similar way.


Here's an example with ChIP-seq data.

A **typical ChIP-seq analysis** of human samples starts from paired-end FASTQ files in the directory `input-dir`:

.. code:: bash

    $ ls /path/to/input-dir/
    my_H3K27ac_sample_R1.fastq.gz  my_H3K27me3_sample_R1.fastq.gz  my_Input_sample_R1.fastq.gz
    my_H3K27ac_sample_R2.fastq.gz  my_H3K27me3_sample_R2.fastq.gz  my_Input_sample_R2.fastq.gz

The :doc:`ChIP-seq` workflow requires the files to be processed via the :doc:`DNA-mapping` workflow first. We therefore run the DNA-mapping workflow :

.. code:: bash

    $ DNA-mapping -i /path/to/input-dir -o /path/to/output-dir --mapq 5 -j 10 --dedup hs37d5

* *--mapq 5* would filter mapped reads for a minimum mapping quality of 5. This would keep only primary alignments from bowtie2, sufficient for downstream analysis.

* *--dedup* would remove PCR duplicates (reads with matching 5' position in the genome), a typical step in ChIP-Seq analysis.

* *-j 10* defines 10 jobs to be run in parallel on the cluster (see below).

* *hs37d5* is the name of the genome (keyword for the yaml). The yaml file corresponding to this genome should exist as `snakePipes/shared/organisms/hs37d5.yaml`. (see :doc:`content/setting_up` for details).

All individual jobs of the workflow will be submitted to the Grid engine using the command specified under /shared/cluster.yaml. The parameter `-j` defines the number of jobs to be run in parallel, while the number of threads per job is hard-coded in the workflows.

**To run the workflow locally**, use the parameter `--local` for local mode and the parameter `-j 10` to specify the maximal number of used CPU threads (here: 10).


Once the DNA-mapping run is finished sucessfully. We can run the ChIP-seq analysis in the same directory. The ChIP-seq workflow would follow up from the DNA-mapping outputs and perform peak calling, create ChIP-input normalized coverage files and also do differential (control-test) analysis if a sample information file is provided (see :ref:`sampleinfo`).

.. code:: bash

    $  ChIP-seq -d /path/to/dna-mapping-output/ hs37d5 chip-samples.yaml



.. _sampleinfo:

The sample_info tsv file
------------------------

Most of the workflows allow users to perform grouped operations as an option, for example
differential expression analysis in RNA-seq workflow, differential binding analysis in
ChIP-Seq workflow, differential open-chromatin analysis in ATAC-seq workflow or merging of
groups in Hi-C workflow. For all this analysis, snakePipes needs a sample_info.tsv file that
contains sample grouping information. In most cases users would want to groups samples by
replicates. The format of the file is as follows:

::

    name   condition
    sample1    control
    sample1    control
    sample2    test
    sample2    test

The name section referes to sample names (without the read suffix), while the condition
section refers to sample group (control/test, male/female, normal/diseased etc..)

Advanced usage
--------------

It's also possible to directly run `snakemake`, please see :doc:`advanced_usage` for more information.
