.. _setting_up:

Setting up snakePipes
=====================

Unlike many other pipelines, setting up snakePipes is easy! All you need is a *linux/OSX system* with *python3-conda* installation.

Installing conda with python3
-----------------------------

Follow the instructions `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`__ to install either
miniconda or anaconda. A minimal version (miniconda) is enough for snakePipes. Get the miniconda installer `here <https://conda.io/miniconda.html>`__.

After installation, check your python path and version :

.. code-block:: bash

    $ which python
    $ /your_path/miniconda3/bin/python

    $ python --version # anything above 3.5 is ok!
    $ Python 3.6.5 :: Anaconda, Inc.

    $ conda --version # only for sanity check
    $ conda 4.5.8

Next, install snakePipes.


Installing snakePipes
---------------------

The easiest way to install snakePipes is via our conda channel. The following command install snakePipes and also creates a conda virtual environment named ``snakePipes``, which you can then activate via ``conda activate snakePipes``.

.. code:: bash

    conda create -n snakePipes -c mpi-ie -c bioconda -c conda-forge snakePipes

This way, the software used within snakePipes do not conflict with the software pre-installed on your terminal or in your python environment.

.. note:: This might take a few minutes depending on the access to conda channels.

Development installation
~~~~~~~~~~~~~~~~~~~~~~~~

If you wish to modify snakePipes you can install it via pip, using our `GitHub repository <https://github.com/maxplanck-ie/snakepipes>`__ .

.. code:: bash

    pip install --user --upgrade git+https://github.com/maxplanck-ie/snakepipes@develop

Instead of providing the URL to ``pip``, you can also `clone <https://help.github.com/articles/cloning-a-repository/>`__ our `GitHub repository <https://github.com/maxplanck-ie/snakepipes>`__ on your computer, and modify the code before running snakePipes. Please see :doc:`advanced_usage` for more information on how to modify and extend snakePipes workflows.

.. note:: There is a difference between installing via conda or installing via pip. The python  installation from user's ``$PATH`` is ignored when installing via conda (first method) while is considered when installing via pip. You must use the ``--develop`` option later when you run ``snakePipes createEnvs``.

.. note:: Using the --user argument would install the program into ``~/.local/bin/``. So make sure to have it in your $PATH before executing any workflow.


Testing whether the installation went fine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After installation, you can activate the snakePipes environment via

.. code:: bash

    conda activate snakePipes

In case you installed conda using the latest version of conda installers (eg. minicoda `4.5.*` or later), the `conda` command might not be available inside an environment. To enable this, export the path to conda/bin in your $PATH (or append the path manually in your `bashrc`)

.. code:: bash

    export PATH="/path/to/miniconda3/bin:$PATH"

Snakemake and pandas are installed along with snakePipes as requirements. Ensure you have them working by testing these commands:

.. code-block:: bash

    snakemake --help
    snakePipes --help


Inspect and modify the setup files
----------------------------------

After installation of snakePipes, all files required to configure it would be installed in a default path. The path to these files can be displayed by running the following command:

.. code:: bash

    snakePipes info

This would show the locations of:

 * **defaults.yaml** Defines default tool and file paths. See :ref:`conda`
 * **cluster.yaml** Defines execution command for the cluster. See :ref:`cluster`
 * **organisms/<organism>.yaml** : Defines genome indices and annotations for various organisms. See :ref:`organisms`
 * Workflow-specific defaults : Defines default options for our command line wrappers. See :ref:`workflowOpts`

You can modify these files to suite your needs before creating the conda environments (see below).


.. _conda:

Install the conda environments
------------------------------

All the tools required for running various pipelines are installed via various conda repositories
(mainly bioconda). The following commands installs the tools and creates the respective conda environments.

.. code:: bash

    snakePipes createEnvs

.. note:: Creating the environments might take 1-2 hours. But it only has to be done once.

.. note::

    ``snakePipes createEnvs`` will also set the ``snakemakeOptions:`` line in the global snakePipes
    ``defaults.yaml`` files. If you have already modified this then use the ``--keepCondaDir`` option.

.. warning::
   If you installed with ``pip`` you must use the ``--develop`` option.

The place where the conda envs are created (and therefore the tools are installed) is defined in ``snakePipes/defaults.yaml``
file on our GitHub repository. You can modify it to suite your needs.

Here are the content of *defaults.yaml*::

    snakemakeOptions: '--use-conda --conda-prefix /data/general/scratch/conda_envs'

.. note::

    Whenever you change the `snakemakeOptions:` line in `defaults.yaml`, you should run
    `snakePipes createEnvs` to ensure that the conda environments are then created.

Running ``snakePipes createEnvs`` is not strictly required, but facilitates multiple users using the same snakePipes installation.


.. _organisms:

Configure the organisms
-----------------------

For each organism of your choice, create a file called ``shared/organisms/<organism>.yaml`` and
fill the paths to the required files next to the corresponding yaml entry. For common organisms, the required files are downloaded and the yaml entries can be created automatically via the workflow ``createIndices``.

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
    # Needed for DNA-mapping workflow
    bowtie2_index: "/data/repository/organisms/dm3_ensembl/BowtieIndex/genome"
    # index of the genome.fasta using HISAT2, needed for RNA-seq workflow
    hisat2_index: "/data/repository/organisms/dm3_ensembl/HISAT2Index/genome"
    # needed by HISAT2 for RNA-seq workflow
    known_splicesites: "/data/repository/organisms/dm3_ensembl/ensembl/release-78/HISAT2/splice_sites.txt"
    bwa_index: "/data/repository/organisms/dm3_ensembl/BWAindex/genome.fa"
    # index of the genome.fasta using STAR, needed for RNA-seq workflow
    star_index: "/data/repository/organisms/dm3_ensembl/STARIndex/"
    # Needed for QC and annotation in DNA-mapping/RNA-Seq workflows
    genes_bed: "/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.bed"
    # Needed for QC and annotation in DNA-mapping/RNA-Seq workflows
    genes_gtf: "/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.gtf"
    # OPTIONAL. For QC and filtering of regions in multiple workflows.
    blacklist_bed:
    # STRING. Name of the chromosomes to ignore for calculation of normalization factors for coverage files
    ignore_forNorm: "U Uextra X XHet YHet dmel_mitochondrion_genome"

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

 - `Human (GRCh38, Gencode release 29) <https://zenodo.org/record/2650763>`__
 - `Mouse (GRCm38/mm10, Gencode release m19) <https://zenodo.org/record/2650854>`__
 - `Mouse (GRCm37/mm9, Gencode release 1) <https://zenodo.org/record/2650849>`__
 - `Fruit fly (dm6, Ensembl release 94) <https://zenodo.org/record/2650762>`__

To use these, simply download and extract them. You will then need to modify the provided YAML file to indicate exactly where the indices are located (i.e., replace ``/data/processing/ryan`` with whatever is appropriate).

.. _cluster:

Configure your cluster
----------------------

The ``cluster.yaml`` file contains both the default memory requirements as well as two options passed to snakemake that control how jobs are submitted to the cluster and files are retrieved::

    snakemake_latency_wait: 300
    snakemake_cluster_cmd: module load slurm; SlurmEasy --mem-per-cpu {cluster.memory} --threads {threads} --log {snakePipes_cluster_logDir} --name {rule}.snakemake 
    snakePipes_cluster_logDir: cluster_logs
    __default__:
        memory: 8G
    snp_split:
        memory: 10G

If you have cloned the repository locally, the file is located under ``snakePipes/shared/``.

You can change the default per-core memory allocation if needed here. Importantly, the ``snakemake_cluster_cmd`` 
option must be changed to match your needs (see table below). Whatever command you specify must include 
a ``{cluster.memory}`` option and a ``{threads}`` option. You can specify other required options here as well. 
The ``snakemake_latency_wait`` value defines how long snakemake should wait for files to appear 
before throwing an error. The default of 300 seconds is typically reasonable when a file system such as 
`NFS <https://en.wikipedia.org/wiki/Network_File_System>`__ is in use. Please also note that there are additional memory 
settings for each workflow in ``snakePipes/workflows/[workflow]/cluster.yaml`` that you might need to adjust. 

``snakePipes_cluster_logDir:`` can be used like a wildcard in `snakemake_cluster_cmd` to specify the directory 
for the stdout and stderr files from a job that is running on the cluster. This is given separate to make sure 
the directory exists before execution. A relative path is treated relative to the ouput directory of the workflow. 
If you want, you can also give an absolute log directory starting with /.

==================== ======================================================================================
 Scheduler/Queuing        snakemake_cluster_cmd example                                                                                                    
==================== ======================================================================================
 **slurm**            .. code:: bash                                                                                       
                                          
                        snakemake_cluster_cmd: module load slurm; sbatch --ntasks-per-node=1 
                           -c {threads} -J {rule}.snakemake --mem-per-cpu={cluster.memory} 
                           -p MYQUEUE -o {snakePipes_cluster_logDir}/{rule}.%j.out 
                           -e {snakePipes_cluster_logDir}/{rule}.%j.err
                        snakePipes_cluster_logDir: cluster_logs
                        
 **PBS/Torque**       .. code:: bash                                                                                       
                                          
                        snakemake_cluster_cmd: qsub -N {rule}.snakemake
                           -q MYQUEUE -l pmem={cluster.memory} 
                           -l walltime=20:00:00 -l nodes=1:ppn={cluster.threads} 
                           -o {snakePipes_cluster_logDir}/{rule}.\$PBS_JOBID.out 
                           -e {snakePipes_cluster_logDir}/{rule}.\$PBS_JOBID.err
                        snakePipes_cluster_logDir: cluster_logs        
                        
 **SGE**              *Please send us a working example!*                
==================== ======================================================================================



.. _workflowOpts:

Configure default options for workflows
---------------------------------------

The default options for all command-line arguments as well as for the cluster (memory) are stored in the workflow-specific folders. If you have cloned the repository locally, these files are located under ``snakePipes/workflows/<workflow_name>`` folder. You can modify the values in these yamls to suite your needs. Most of the default values could also be replaced from the command line wrappers while executing a workflow.

Below are some of the workflow defaults from the DNA-mapping pipeline. Empty sections means no default is set:

.. parsed-literal::
    ## key for the genome name (eg. dm3)
    genome:
    ## FASTQ file extension (default: ".fastq.gz")
    ext: '.fastq.gz'
    ## paired-end read name extension (default: ['_R1', "_R2"])
    reads: [_R1, _R2]
    ## mapping mode
    mode: mapping
    mapping_prg: Bowtie2
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
    gcbias: false
    ## Retain only de-duplicated reads/read pairs
    dedup: false
    ## Retain only reads with at least the given mapping quality
    mapq: 0

Test data
---------

Test data for the various workflows is available at the following locations:

 - `DNA mapping <https://zenodo.org/record/1346303>`__
 - `ChIP-seq <https://zenodo.org/record/2624281>`__
 - `ATAC-seq <https://zenodo.org/record/2624323>`__
 - `RNA-seq <https://zenodo.org/record/2624408>`__
 - `HiC <https://zenodo.org/record/2624479>`__
 - `WGBS <https://zenodo.org/record/2624498>`__
 - `scRNA-seq <https://zenodo.org/record/2624518>`__
