[![linux](https://github.com/maxplanck-ie/snakepipes/actions/workflows/linux.yml/badge.svg)](https://github.com/maxplanck-ie/snakepipes/actions/workflows/linux.yml)
[![osx](https://github.com/maxplanck-ie/snakepipes/actions/workflows/osx.yml/badge.svg)](https://github.com/maxplanck-ie/snakepipes/actions/workflows/osx.yml)
[![pytest](https://github.com/maxplanck-ie/snakepipes/actions/workflows/pytest.yml/badge.svg)](https://github.com/maxplanck-ie/snakepipes/actions/workflows/pytest.yml)
[![readthedocs](https://readthedocs.org/projects/snakepipes/badge/?version=latest)](https://snakepipes.readthedocs.io/en/latest/)
[![citation](https://zenodo.org/badge/54579435.svg)](https://zenodo.org/badge/latestdoi/54579435)

# SnakePipes

snakePipes are flexible and powerful workflows built using [snakemake](https://github.com/snakemake/snakemake) that simplify the analysis of NGS data.
![snakePipes](.docs/content/images/snakePipes_small.png)

## Workflows

- DNAmapping*  
- ChIPseq*  
- mRNAseq*  
- ncRNAseq*  
- ATACseq*  
- scRNAseq  
- HiC  
- makePairs*  
- Whole Genome Bisulfite Seq/WGBS  

(* also available in allele-specific mode)

## Installation

[Conda](https://docs.conda.io/en/latest/#) is a pre-requisite for snakePipes. So make sure this is [installed](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) before.

Afterwards you can create a snakePipes environment containing the installation by running:

 > conda create -n snakepipes -c mpi-ie -c bioconda -c conda-forge snakePipes

In case you'd like a development version, you can install snakePipes directly from github using pip:

 > git clone git@github.com:maxplanck-ie/snakepipes.git  
 > cd snakepipes
 > pip install .

Make sure the environment you are installing this version into has python version 3.11 or later.

After the installation some configurations have to be set, for which we refer to the documentation.

## Documentation

For detailed documentation on setup and usage, please visit the [documentation](https://snakepipes.readthedocs.io/en/latest/).

## Citation

If you adopt/run snakePipes for your analysis, please cite it as follows :

Bhardwaj, Vivek, Steffen Heyne, Katarzyna Sikora, Leily Rabbani, Michael Rauer, Fabian Kilpert, Andreas S. Richter, Devon P. Ryan, and Thomas Manke. 2019. “snakePipes: Facilitating Flexible, Scalable and Integrative Epigenomic Analysis.” Bioinformatics , May. [doi:10.1093/bioinformatics/btz436](https://doi.org/10.1093/bioinformatics/btz436).

## Note

SnakePipes are under active development. We appreciate your help in improving it further. Please use issues to the GitHub repository for feature requests or bug reports.
