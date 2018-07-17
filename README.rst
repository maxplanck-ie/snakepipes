===========================================================
MPI-IE Snakemake workflows : snakePipes
===========================================================

.. image:: https://readthedocs.org/projects/snakepipes/badge/?version=latest
    :target: http://snakepipes.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.org/maxplanck-ie/snakepipes.svg?branch=develop
    :target: https://travis-ci.org/maxplanck-ie/snakepipes
    :alt: Build Staus

.. image:: https://zenodo.org/badge/54579435.svg
    :target: https://zenodo.org/badge/latestdoi/54579435
    :alt: Citation

snakePipes are our flexible and powerful workflows built using `snakemake <snakemake.readthedocs.io>`__ that simplify the analysis of NGS data.

Workflows available
--------------------

- DNA-mapping (normal and allele-specific)
- ChIP-seq (normal and allele-specific)
- RNA-seq (normal and allele-specific)
- ATAC-seq (normal and allele-specific)
- scRNA-seq
- Hi-C
- Whole Genome Bisulfite Seq/WGBS

Installation
-------------

Snakepipes uses conda for installation and dependency resolution, so you will need to `install conda <https://conda.io/docs/user-guide/install/index.html>`__ first.

Afterward, simply run the following:

    conda create -n snakePipes -c mpi-ie -c bioconda -c conda-forge snakePipes

This will create a new conda environment called "snakePipes" into which snakePipes is installed. You will then need to create the conda environments needed by the various workflows. To facilitate this we provide the `snakePipes` command:

1. `source activate snakePipes` to activate the appropriate conda environment.
2. `snakePipes createEnvs` to create the various environments and register GATK.

A few additional steps you can then take:

1. **Modify/remove/add the organism yaml files appropriately** : these yaml files would contain location of appropriate
GTF files and genome indexes corresponding to different organisms. The location of these files after installation can be
found using `snakePipes --info` command.

2. **Modify the cluster.yaml file appropriately** : This yaml file contains information about your cluster scheduler (SGE/slurm).
Location revealed using `snakePipes --info` command.


Documentation
--------------

For detailed documentation on setup and usage, please visit our `read the docs page <https://snakepipes.readthedocs.io/en/latest/>`__.


Citation
-------------

If you adopt/run snakePipes for your analysis, cite it as follows :

    snakePipes: flexible and scalable NGS analysis pipelines built using snakemake. The MPI-IE Bioinformatics Facility. http://doi.org/10.5281/zenodo.1146540


Notice
-------------

SnakePipes are currently unstable and under active development.
