Setting up snakePipes
==========================

Unlike many other pipelines, setting up snakePipes is easy! All you need is a *conda* installation with *python3*.

Install conda with python3
---------------------------

.. note:: If you have any other, non-conda version/installation of python configured in your $PATH,
        ~/.bashrc or ~/.bash_profile, please remove them first in order to avoid conflicts.

Then follow the instructions [here](https://conda.io/docs/user-guide/install/index.html) to install either
miniconda or anaconda. A minimal version (miniconda) is enough for snakePipes. Get miniconda installer [here](https://conda.io/miniconda.html).

After installation, check your python path and version :

.. code:: bash

    $ which python
    $ /your_path/miniconda3/bin/python

    $ python --version # anything above 3.6 is ok!
    $ Python 3.6.5 :: Anaconda, Inc.

    $ conda --version # only for sanity check
    $ conda 4.5.8

Next, install snakePipes.


Install snakePipes
-------------------------

The easiest way to install snakePipes is via our conda channel. The following command also creates a
conda virtual environment named `snakePipes`, which you can then activate via `source activate snakePipes`.

.. code:: bash
    conda create -n snakePipes -c mpi-ie -c bioconda -c conda-forge snakePipes

Another way is via pip, using our [GitHub repository](https://github.com/maxplanck-ie/snakepipes).

.. code:: bash
    pip install --user --upgrade git+https://github.com/maxplanck-ie/snakepipes@develop

Finally, for advanced users who want to modify the source code, you can install via pip,
after cloning the source code from our [GitHub repository](https://github.com/maxplanck-ie/snakepipes).

.. code:: bash
    git clone https://github.com/maxplanck-ie/snakepipes.git && cd snakePipes
    # then configure as per your wish
    # now install
    pip install --user --upgrade /path/to/local/snakepipes/

.. note:: There is a difference between installing via conda or installing via pip. The python installation from user's
$PATH is ignored when installing via conda (first method) while is considered when installing via pip.

.. note:: Using the --user argument would install the program into `~/.local/bin/`. So make sure to have it in your $PATH

.. code:: bash
    export PATH=~/.local/bin:$PATH

Snakemake and pandas are installed as requirements to snakePipes. Ensure you have everything working by testing these commands:

.. code:: bash
    snakemake --help
    snakePipes --help


Install the conda environments
--------------------------------

All the tools required for running various pipelines are installed via various conda repositories
(mainly bioconda). The following commands installs the tools and creates the respective conda environments.

.. code:: bash
    snakePipes createEnvs

.. note:: Creating the environments might take 1-2 hours. But it only has to be done once.

The place where the conda envs are created (and therefore the tools are installed) is defined in `snakePipes/defaults.yaml`
file on our GitHub repository. You can modify it to suite your needs.

Here's the content of *defaults.yaml*:

.. parsed-literal::
    snakemake_options: '--use-conda --conda-prefix /data/general/scratch/conda_envs'
    tempdir: /data/extended/

The `tempdir` path could be changed to any suitable directory that can hold the temporary files during pipeline execution.


Configure the organisms
----------------------------

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

Not all files are required for all pipelines, but we recommend to keep all required files ready nevertheless.


################################################################################
# Call snakemake directly, i.e. without using the wrapper script:
#
# Please save a copy of this config yaml file and provide an adjusted config
# via '--configfile' parameter!
# example call:
#
# snakemake --snakefile /path/to/snakemake_workflows/workflows/ATAC-seq/Snakefile
#           --configfile /path/to/snakemake_workflows/workflows/ATAC-seq/defaults.yaml
#           --directory /path/to/outputdir
#           --cores 32
################################################################################
