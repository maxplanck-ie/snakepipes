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
The initial dependency for snakePipes is conda with python3. Other dependencies (snakemake and pandas) are installed during installation of snakePipes.
Install snakePipes from our GitHub repo as follows:

    pip install --user --upgrade git+https://github.com/maxplanck-ie/snakepipes@master

snakePipes is now installed, however you will need to perform a few steps to configure snakePipes and finish the setup:

1. *Modify/remove/add the organism yaml files appropriately* : these yaml files would contain location of appropriate
GTF files and genome indexes corresponding to different organisms. The location of these files after installation can be
found using `snakePipes --info` command.

2. *Modify the cluster.yaml file appropriately* : This yaml file contains information about your cluster scheduler (SGE/slurm).
Location revealed using `snakePipes --info` command.

3. *Create the conda environments* : conda environments for workflows could be created using `snakePipes --createEnvs` command.
This takes a while, but only need to be done once. Optionally, workflow-specific env could also be created automatically when you 
run any specific workflow for the first time.


.. note:: `snakePipes createEnvs` will also set the `snakemake_options:` line in the global snakePipes `defaults.yaml` files. If you have already modified this then use the `--keepCondaDir` option.

.. note:: Whenever you change the `snakemake_options:` line in `defaults.yaml`, you should run `snakePipes createEnvs` to ensure that the conda environments are then created.

.. note:: Running `snakePipes createEnvs` is not strictly required, but facilitates multiple users using the same snakePipes installation.

Documentation
--------------

For documentation on setup and usage, please visit our `read the docs page <https://snakepipes.readthedocs.io/en/latest/>`__.

Citation
-------------

If you adopt/run snakePipes for your analysis, cite it as follows :

```

snakePipes: flexible and scalable NGS analysis pipelines built using snakemake. The MPI-IE Bioinformatics Facility. http://doi.org/10.5281/zenodo.1146540

```

Notice
-------------

SnakePipes are currently unstable and under active development.
