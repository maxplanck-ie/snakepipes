===========================================================
MPI-IE Snakemake workflows : snakePipes
===========================================================

.. image:: https://readthedocs.org/projects/snakepipes/badge/?version=latest
    :target: http://snakepipes.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.org/maxplanck-ie/snakepipes.svg?branch=develop
    :target: https://travis-ci.org/maxplanck-ie/snakepipes
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
- ATAC-seq (normal and allele-specific)
- scRNA-seq
- Hi-C

Installation
-------------

Currently the workflows don't require installation of the scripts (install.sh simply creates the link to available workflows).
Simply clone this github repository and modify the configuration files under `shared/organisms` and `shared/paths.yaml` with your own paths.

### Points to consider while adopting snakePipes for your institution :
 - While some paths under `shared/paths.yaml` are loaded directly using abspath, other are loaded via `module load`. This is how we have
   implemented this so far for our case. Replace these paths with your own method of loading programs (eg. via source activate).
 - Python 3 and Snakemake are the absolute dependencies. Snakemake itself is loaded via module load in our case. Replace and install as per your
   configuration.

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
