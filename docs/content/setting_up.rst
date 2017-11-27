Setting up the workflows
==========================

To setup snakePipes after a fresh download, you need to do the following.

Set up slurm and snakemake
--------------------------

The pipelines require snakemake in order to work, and slurm in order to submit jobs to the cluster.
If you don't have slurm configured with the cluster, you can skip this and run the pipelines locally using the
`--local` option in the wrappers.

Edit the paths to the required programs
---------------------------------------

The paths to the required programs can be configured under `shared/paths.yaml`. This contains a list of all
programs required, but not all workflows required all of these programs to be installed and therefore some of
them can be skipped depending upon the workflow used.

.. warning:: Do not edit the yaml keywords corresponding to each required entry.

Configure the organisms
------------------------

For each organism of your choice, create a file called `shared/organisms/<organism>.yaml` and
fill the paths to the required files next to the corresponding yaml entry.

.. warning:: Do not edit the yaml keywords corresponding to each required entry.

An example from drosophila genome dm3 is below.

.. parsed-literal::

    genome_size: 142573017
    genome_fasta: "/data/repository/organisms/dm3_ensembl/genome_fasta/genome.fa"
    genome_index: "/data/repository/organisms/dm3_ensembl/genome_fasta/genome.fa.fai"
    genome_2bit: "/data/repository/organisms/dm3_ensembl/genome_fasta/genome.2bit"
    bowtie2_index: "/data/repository/organisms/dm3_ensembl/BowtieIndex/genome"
    hisat2_index: "/data/repository/organisms/dm3_ensembl/HISAT2Index/genome"
    bwa_index: "/data/repository/organisms/dm3_ensembl/BWAindex/genome.fa"
    known_splicesites: "/data/repository/organisms/dm3_ensembl/ensembl/release-78/HISAT2/splice_sites.txt"
    star_index: "/data/repository/organisms/dm3_ensembl/STARIndex/"
    genes_bed: "/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.bed"
    genes_gtf: "/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.gtf"
    blacklist_bed:
    ignore_forNorm: "U Uextra X XHet YHet dmel_mitochondrion_genome"

Not all files are required for all workflows, but we recommend to keep all required files ready nevertheless..
