.. _setting_up:

Setting up snakePipes
=====================

Unlike many other pipelines, setting up snakePipes is easy! All you need is a *linux/OSX system* with a working *conda* installation.

Installing conda
----------------

.. note::
   Latest snakePipes versions (3.0 onwards) require conda >= 23.10.0 If you have an older version of conda, please don't use it. You may try installing the extension (conda-libmamba-solver) in your base environment, and setting this as default (conda config --set solver libmamba) but this has not been tested by us, hence YMMV. Mamba used to be a pre-requisite for snakePipes (versions 2.5.3 through 2.8.1), refer to older docs in such versions.

Follow the instructions `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`__ to install either miniconda or anaconda first.
After installation, check your python path and version :

.. code-block:: bash

    $ command -v python
    $ <your_chosen_installation_path>/bin/python

    $ python --version # anything above 3.5 is ok!
    $ Python 3.6.5 :: Anaconda, Inc.

Installing snakePipes
---------------------

The easiest way to install snakePipes is via our conda channel. The following command install snakePipes and also creates a conda virtual environment named ``snakePipes``, which you can then activate via ``conda activate snakePipes``. Specifying snakePipes version avoids issues with conda's environment solver.

.. code:: bash

    conda create -n snakePipes -c mpi-ie -c conda-forge -c bioconda snakePipes

This way, the software used within snakePipes do not conflict with the software pre-installed on your terminal or in your python environment.

Now, we should activate this environment:

.. code:: bash

    conda activate snakePipes


Alternatively, snakePipes can also be installed using pip. You can clone a branch (given that it's version 3.0.0 or later) and just install with pip:

.. code:: bash

    pip install .

Just make sure you have python 3.11 or later (cap by snakemake) in your environment.
In case you'd like to develop snakePipes, extra dependencies for the documentation:

.. code:: bash

    pip install .[docs]

or for the pytests:

.. code:: bash

    pip install .[actions]


Configuring snakePipes
----------------------

Finally, config files ``defaults.yaml`` and snakemake profile should be modified to match your compute infrastructure. The location of this file can be found out by executing:

.. code:: bash

    snakePipes info

This would return you where the global configuration file is located.
Two fields are important to set:

++++++++++++++++
snakemakeProfile
++++++++++++++++
Defines a `snakemake profile <https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`__ to use.
By default this translates to a pre-shipped 'local' profile (and points to a location relative to the snakePipes package directory).
The local profile runs all jobs without a submission system. 

Another profile shipped within the repository is a default slurm profile (using snakemake-executor-plugin-cluster-generic). 
In case you want to use this you can set the snakemakeProfile value to ``shared/profiles/snakepipes_genericprofile``.
After changing the value of snakemakeProfile, you should re-run ``snakePipes info``, which will also print out the full directory of the profile used.
If you want to use the snakepipes_genericprofile, make sure to review the following entries in the profile yaml file with respect to your infrastructure:

 * ``module load slurm &&`` - could be omitted
 * ``resources.partition`` - set to your slurm partition
 * ``conda-prefix`` - set to your preferred location where snakePipes environments should be stored. You can set this value by running `snakePipes config --condaEnvDir` and providing the respective path.
 * ``resources`` - make sure default resources make sense for your infrastructure
 * ``ccancel.sh`` - refers to the ccancel.sh file inside the profile directory and contains instructions on how to kill submitted jobs (on failure / interruption of snakemake). The module command could be omitted here as before

 In case you are using your own snakemake profile already, you can define them here as well. Acceptable values in snakemakeProfile are:
 
 * absolute path to a snakemake profile directory
 * a relative path to a snakemake profile (relative to the package directory)
 * The name of a `global snakemake profile <https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`__ 

 If you use your own profile, just make sure that at least these values are set in your profile:

 * use-conda: true
 * conda-prefix: /path/to/prefix
 * conda-frontend: conda

Additionaly, rule resources are defined in the pre-shipped profiles. 
In case you use your own you'd want to have these set in your profile as well.

+++++++
tempDir
+++++++
The temp directory to use. Defaults to /scratch/local.

After setting the defaults, the conda environments can be created. 

.. _conda:

Create the conda environments
-----------------------------

All the tools required for running various pipelines are installed via various conda repositories
(mainly bioconda). The following commands installs the tools and creates the respective conda environments.
Note that the conda-prefix is defined in your profile (and defaults to /tmp). Thus make sure you have set your profile appropriately.
It is important that the conda-prefix is a location that is accessible by your compute nodes as well.
Finally, make sure you have a conda installation with libmamba as the solver (conda version 23.10.0 or later), as this speeds up the process.

.. code:: bash

    snakePipes createEnvs


.. _organisms:

Configure the organisms
-----------------------

For each organism of your choice, create a file called ``<organism>.yaml`` in the folder specified by ``organismsDir`` in **defaults.yaml** and
fill the paths to the required files next to the corresponding yaml entry. For common organisms, the required files are downloaded and the yaml entries can be created automatically via the workflow ``createIndices``.

Note that the organism yamls that come with the installation are only appropriate internally for MPI-IE, and as an external you need to create / download your own.

The yaml files look like this after the setup (an example from drosophila genome ``dm3``) :

.. parsed-literal::

    # Integer, size of genome in base-pairs
    genome_size: 142573017
    # path to genome.fasta for mapping
    genome_fasta: "/data/repository/organisms/dm3_ensembl/genome_fasta/genome.fa"
    # path to genome.fasta.fai (fasta index) for mapping
    genome_index: "/data/repository/organisms/dm3_ensembl/genome_fasta/genome.fa.fai"
    # OPTIONAL. Needed for GC bias estimation by deepTools
    genome_2bit: "/data/repository/organisms/dm3_ensembl/genome_fasta/genome.2bit"
    # Needed for DNAmapping workflow
    bowtie2_index: "/data/repository/organisms/dm3_ensembl/BowtieIndex/genome"
    # index of the genome.fasta using HISAT2, needed for RNA-seq workflow
    hisat2_index: "/data/repository/organisms/dm3_ensembl/HISAT2Index/genome"
    # needed by HISAT2 for RNA-seq workflow
    known_splicesites: "/data/repository/organisms/dm3_ensembl/ensembl/release-78/HISAT2/splice_sites.txt"
    bwa_index: "/data/repository/organisms/dm3_ensembl/BWAindex/genome.fa"
    # index of the genome.fasta using STAR, needed for RNA-seq workflow
    star_index: "/data/repository/organisms/dm3_ensembl/STARIndex/"
    # Needed for QC and annotation in DNAmapping/RNA-Seq workflows
    genes_bed: "/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.bed"
    # Needed for QC and annotation in DNAmapping/RNA-Seq workflows
    genes_gtf: "/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.gtf"
    # OPTIONAL. For QC and filtering of regions in multiple workflows.
    blacklist_bed:
    # STRING. Name of the chromosomes to ignore for calculation of normalization factors for coverage files
    ignoreForNormalization: "U Uextra X XHet YHet dmel_mitochondrion_genome"

.. warning:: Do not edit the yaml keywords corresponding to each required entry.

.. note:: Some fields are optional and can be left empty. For example, if a blacklist file
          is not available for your organism of interest, leave `blacklist_bed:` empty.
          Files for either STAR or HISAT2 could be skipped for RNA-seq if the respective
          aligner is not used. We nevertheless recommended providing all the files, to allow
          more flexible analysis.

After setting up the yamls, we can execute a snakePipes workflow on the organism of choice by referring to the **organism** as ``dm3``, where the keyword **dm3** matches the name of the yaml file (dm3.yaml).

.. note:: The name of the yaml file (except the .yaml suffix) is used as keyword to refer to the organism while running the workflows.

Download premade indices
------------------------

For the sake of convenience, we provide premade indices for the following organisms:

 - `Human (GRCh38, Gencode release 29) <https://zenodo.org/record/4471116>`__
 - `Mouse (GRCm38/mm10, Gencode release m19) <https://zenodo.org/record/4468065>`__
 - `Mouse (GRCm37/mm9, Gencode release 1) <https://zenodo.org/record/4478284>`__
 - `Fruit fly (dm6, Ensembl release 94) <https://zenodo.org/record/4478414>`__

To use these, simply download and extract them. You will then need to modify the provided YAML file to indicate exactly where the indices are located (i.e., replace ``/data/processing/ryan`` with whatever is appropriate).




.. _workflowOpts:

Configure default options for workflows
---------------------------------------

The default options for all command-line arguments as well as for the cluster (memory) are stored in the workflow-specific folders. If you have cloned the repository locally, these files are located under ``snakePipes/workflows/<workflow_name>`` folder. You can modify the values in these yamls to suite your needs. Most of the default values could also be replaced from the command line wrappers while executing a workflow.

Below are some of the workflow defaults from the DNAmapping pipeline. Empty sections means no default is set:

.. parsed-literal::
    ## key for the genome name (eg. dm3)
    genome:
    ## FASTQ file extension (default: ".fastq.gz")
    ext: '.fastq.gz'
    ## paired-end read name extension (default: ['_R1', "_R2"])
    reads: [_R1, _R2]
    ## mapping mode
    mode: mapping
    aligner: Bowtie2
    ## Number of reads to downsample from each FASTQ file
    downsample:
    ## Options for trimming
    trim: False
    trimmer: cutadapt
    trimmerOptions:
    ## Bin size of output files in bigWig format
    bwBinSize: 25
    ## Run FASTQC read quality control
    fastqc: false
    ## Run computeGCBias quality control
    GCBias: false
    ## Retain only de-duplicated reads/read pairs
    dedup: false
    ## Retain only reads with at least the given mapping quality
    mapq: 0

Test data
---------

Test data for the various workflows is available at the following locations:

 - `DNAmapping <https://zenodo.org/record/3707259>`__
 - `ChIPseq <https://zenodo.org/record/2624281>`__
 - `ATACseq <https://zenodo.org/record/3707666>`__
 - `mRNAseq <https://zenodo.org/record/3707602>`__
 - `ncRNAseq <https://zenodo.org/deposit/3707749>`__
 - `HiC <https://zenodo.org/record/3707714>`__
 - `WGBS <https://zenodo.org/record/3707727>`__
 - `scRNAseq <https://zenodo.org/record/3707747>`__
