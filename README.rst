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

snakePipes are our in-house, flexible and powerful workflows built using `snakemake <snakemake.readthedocs.io>`__ that simplify the analysis of NGS data.

Workflows available
--------------------

- DNA-mapping (normal and allele-specific)
- ChIP-seq (normal and allele-specific)
- RNA-seq (normal and allele-specific)
- ATAC-seq (normal and allele-specific)
- scRNA-seq
- Hi-C

Installation
-------------

Snakepipes uses conda for installation and dependency resolution, so you will need to `install conda <https://conda.io/docs/user-guide/install/index.html>`__ first. Then create an environment for snakepipes:

    conda create -n snakepipes-1.0.0 -c mpi-ie -c bioconda -c conda-forge snakepipes==1.0.0

This snakepipes is then installed, however you will need to perform a few steps to actually configure snakepipes and finish the setup:

  1. Modify/remove/add the organism yaml files appropriately
  2. Modify the cluster.yaml file appropriately
  3. Create the per-workflow conda environments

To facilitate this, you can use the `snakePipes info` command to list the locations of the various yaml files. Edit or remove these to suite your needs and then run `snakePipes createEnvs` to create the per-workflow conda environments. Note that this will also set the `snakemake_options:` line in the various workflow `default.yaml` files, by default. If you have already modified this then use the `--keepCondaDir` option.

.. note:: Whenever you change the `snakemake_options:` line in a workflow `default.yaml`, you should run `snakePipes createEnvs` to ensure that the conda environments are then created.

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
