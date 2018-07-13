Running snakePipes
=======================

Pipelines under snakePipes are designed in a modular way such that all are configured and ran in a similar way.


Here's an example with ChIP-seq data.

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

Here, hs37d5 is the name of the genome. The yaml file corresponding to this genome should exist under `snakePipes/shared/organisms/hs37d5.yaml`.
This yaml file should have paths to the required genome fasta, index, GTF and other annotations (see **Genome configuration file** below).

All individual jobs of the workflow will be submitted to the Grid engine using the command specified under /shared/cluster.yaml.
To run the workflow locally, use the parameter `--local` for local mode and the parameter `-j 48` to specify the maximal
number of used CPU threads (here: 48) or concurrent running Slurm jobs (actual used threads are defined in each rule).


Genome configuration file
----------------------------

Setting up of **organism yaml** files has been described in :doc:`/content/setting_up.rst`.
Below is an explanation of each key mentioned in the organism yaml file.

    genome_size: # Integer, size of genome in base-pairs
    genome_fasta: # path to genome.fasta for mapping
    genome_index: # path to genome.fasta.fai (fasta index) for mapping
    genome_2bit: # OPTIONAL. Needed for GC bias estimation by deepTools
    bowtie2_index: # Needed for DNA-mapping workflow
    hisat2_index: # index of the genome.fasta using HISAT2, needed for RNA-seq workflow
    known_splicesites: # needed by HISAT2 for RNA-seq workflow
    star_index: # index of the genome.fasta using STAR, needed for RNA-seq workflow
    genes_bed: # Needed for QC and annotation in DNA-mapping/RNA-Seq workflows
    genes_gtf: # Needed for QC and annotation in DNA-mapping/RNA-Seq workflows
    blacklist_bed: # OPTIONAL. For QC and filtering of regions in multiple workflows.


.. note:: Some fields are optional and can be left empty. For example, if a blacklist file
          is not available for your organism of interest, leave `blacklist_bed:` empty.
          Files for either STAR or HISAT2 could be skipped for RNA-seq if the respective
          aligner is not used. We nevertheless recommended providing all the files, to allow
          more flexible analysis.
