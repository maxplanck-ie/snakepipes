===========================================================
snakePipes
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


snakePipes are flexible and powerful workflows built using `Snakemake <https://snakemake.readthedocs.io>`__ that simplify the analysis of NGS data.

.. image:: ./docs/content/images/snakePipes_small.png
   :scale: 20 %
   :width: 100 px
   :height: 100 px
   :align: right

Workflows available
--------------------

- DNAmapping*
- ChIPseq*
- mRNAseq*
- ncRNAseq*
- ATACseq*
- scRNAseq
- HiC
- Whole Genome Bisulfite Seq/WGBS

**(*Also available in "allele-specific" mode)**

Installation
-------------

Snakepipes is a set of Snakemake workflows which use conda for installation and dependency resolution, so you will need to `install conda <https://conda.io/docs/user-guide/install/index.html>`__ first.

Afterward, simply run the following:

``conda install mamba -c conda-forge && mamba create -n snakePipes -c mpi-ie -c bioconda -c conda-forge snakePipes``

This will create a new conda environment called "snakePipes" into which snakePipes is installed. You will then need to create the conda environments needed by the various workflows. To facilitate this we provide the `snakePipes` commands:

* ``conda activate snakePipes`` to activate the appropriate conda environment.
* ``snakePipes createEnvs`` to create the various environments.

Indices and annotations needed to run the workflows could be created by a simple command :

``createIndices --genomeURL <path/url to your genome fasta> --gtfURL <path/url to genes.gtf> -o <output_dir> <name>``

where `name` refers to the name/id of your genome (specify as you wish).

A few additional steps you can then take:

1. **Modify/remove/add the organism yaml files appropriately** : these yaml files would contain location of appropriate
GTF files and genome indexes corresponding to different organisms. The location of these files after installation can be
found using ``snakePipes info`` command.

2. **Modify the cluster.yaml file appropriately** : This yaml file contains settings for your cluster scheduler (SGE/slurm).
Location revealed using ``snakePipes info`` command.


Documentation
--------------

For detailed documentation on setup and usage, please visit our `read the docs page <https://snakepipes.readthedocs.io/en/latest/>`__.


Citation
-------------

If you adopt/run snakePipes for your analysis, cite it as follows :

Bhardwaj, Vivek, Steffen Heyne, Katarzyna Sikora, Leily Rabbani, Michael Rauer, Fabian Kilpert, Andreas S. Richter, Devon P. Ryan, and Thomas Manke. 2019. “snakePipes: Facilitating Flexible, Scalable and Integrative Epigenomic Analysis.” Bioinformatics , May. `doi:10.1093/bioinformatics/btz436 <https://doi.org/10.1093/bioinformatics/btz436>`__


Note
-------------

SnakePipes are under active development. We appreciate your help in improving it further. Please use issues to the GitHub repository for feature requests or bug reports.
