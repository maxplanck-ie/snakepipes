===========================================================
MPI-IE Snakemake workflows : snakePipes
===========================================================

.. image:: https://readthedocs.org/projects/snakepipes/badge/?version=latest
    :target: http://snakepipes.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.org/maxplanck-ie/snakemake_workflows.svg?branch=develop
    :target: https://travis-ci.org/maxplanck-ie/snakemake_workflows
    :alt: Build Status

.. image:: https://zenodo.org/badge/54579435.svg
    :target: https://zenodo.org/badge/latestdoi/54579435
    :alt: Citation

snakePipes are our in-house, flexible and powerful workflows built using `snakemake <snakemake.readthedocs.io>`__ that simplify the analysis of NGS data.

Workflows available
--------------------

- DNA-mapping (normal and allele-specific)
- ChIP-seq (normal and allele-specific)
- RNA-seq (normal and allele-specific)
- scRNA-seq
- Hi-C

Installation
-------------

Currently the workflows don't require installation of the scripts. Simply clone this github repository
and modify the configuration files under `shared/organisms` and `shared/paths.yaml` with your own paths.

Documentation
--------------

For documentation on setup and usage, please visit our `read the docs page <https://snakepipes.readthedocs.io/en/latest/>`__.
