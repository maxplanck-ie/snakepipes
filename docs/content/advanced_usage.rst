Advanced usage of snakePipes
============================

snakePipes is designed in a modular way, such that it's easy to modify or extend it. Advanced users or developers can either use the underlying snakemake and Rscripts directly, or extend/add modules to the existing workflows.

Understanding snakePipes implementation
----------------------------------------

The implementation of snakePipes modules has been described in `our preprint <https://www.biorxiv.org/content/early/2018/09/18/407312>`__. Please clone our github repository locally to understand this organisation. Since snakePipes is dependent on `snakemake <https://snakemake.readthedocs.io/en/stable/>`__ and `conda <https://conda.io>`__, we recommend being familier with them first, by reading the documentation. Also, we utilize `bioconda <https://www.nature.com/articles/s41592-018-0046-7>`__ as a source of our biology-related tools implemented in snakePipes.

* `Getting started with snakemake <https://slides.com/johanneskoester/snakemake-short#/>`__
* `Getting started with conda <https://conda.io/docs/user-guide/overview.html>`__
* `How conda is used with snakemake <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html>`__
* `What is bioconda <https://bioconda.github.io/>`__

Once you are familier with snakemake and conda/bioconda, we can look at how snakePipes workflows are implemented.

snakePipes folders
~~~~~~~~~~~~~~~~~~~~~~~~~~

All files needed to be modified in order to extend/modify a workflow, are available under the **snakePipes** directory (``snakepipes/snakePipes``). Here is the structure of this directory::

    .
    ├── common_functions.py
    ├── __init__.py
    ├── parserCommon.py
    ├── shared
    │   ├── cluster.yaml
    │   ├── defaults.yaml
    │   ├── organisms
    │   ├── rscripts
    │   ├── rules
    │   └── tools
    └── workflows
        ├── ATACseq
        ├── ChIPseq
        ├── createIndices
        ├── DNAmapping
        ├── HiC
        ├── mRNAseq
        ├── ncRNAseq
        ├── preprocessing
        ├── scRNAseq
        └──WGBS

* **common_functions.py** contains functions that directly operate on the variables received via various wrappers.
* **parserCommon.py** contains common command-line arguments for the wrappers.


* **shared**: This folder contains some important files.

    * **cluster.yaml**: defines the command for the execution of rules on a cluster or cloud, this command is passed on to the call to snakemake via the wrappers.

    * **defaults.yaml**: defined the default options for snakemake and also defines the temporary directory to store intermediate files.

    * **organisms**: This folder contains yaml files corresponding to each organism (see :doc:`setting_up` for details)

    * **rscripts**: Contains the R wrappers that are invoked via the rules. You would find the Rscripts for DESeq, CSAW and other R pakages here.

    * **rules**: These are the snakemake rules which are invoked during execution of a workflow. Depending upon the global parameters passed on from the wrappers, a rule may/may not be included in a workflow (controlled by various ``if`` conditionals).

    .. warning:: Some rules are shared via multiple workflows, therefore be sure to check each Snakemake file for each workflow to see which rules you need to modify.

    * **tools**: This folder contains online tools which can not be distributed via bioconda, and therefore are included with snakePipes package itself.

* **workflows**: This folder contains files which are specific to each workflow. Under each folder (named by the workflow), you would find a common set of files.

    * **<workflow_name>.py**: The command line python **wrappers** that are visible to users after installation.

    * **Snakefile**: This is the file that collects outputs from various rules, therefore contains the ``rule all`` for each workflow. This Snakefile also controls which rules from ``shared/rules`` folder are included in the final workflow, depending on the global parameters passed on from the wrappers.

    * **internals.snakefile**: Contains some python functions which are specific to the workflow, and therefore can't be included under ``common_functions.py``, these functions are imported in the ``Snakefile``

    * **cluster.yaml** and **defaults.yaml** : contains workflow-specific options for the cluster, and for the wrappers. Modify them to suite your needs.


Calling snakemake directly using the snakefiles
------------------------------------------------

It's possible to directly run ``snakemake`` using the ``Snakefile`` provided in each workflow, therefore surpassing the command-line wrappers. 
In order to do that, you can begin with a copy of ``<workflow_name>.defaults.yaml`` file that you will find in your output folder after running the workflow with --snakemakeOptions ' --dryrun ' and add or adjust further options in that file. This file will contained a merged dictionary from the workflow defaults as well as from the global (cross-workflows) defaults.

Finally, provide an adjusted config via ``--configfile`` parameter to snakemake!

example call::

    snakemake --snakefile /path/to/snakemake_workflows/workflows/ATACseq/Snakefile
              --configfile /path/to/(snakemake_workflows/workflows/ATACseq/)defaults.yaml
              --directory /path/to/outputdir
              --cores 32


Executing the Rscript wrappers outside snakePipes
--------------------------------------------------

It's also possible to use one of our Rscript wrappers present under the ``shared/rscripts`` folder. In order to do that, check how the parameters are supplied to the wrappers in the corresponding rule.

For example, in order to execute the DESeq2 wrapper, we can look at how it's done via the DESeq2 rule under ``shared/rules/DESeq2.Snakefile``

example call::

    cd DeSeq2_test &&
    Rscript /path/to/shared/rscripts/DESeq2.R \
    ${input.sample_info} \
    ${input.counts_table} \
    ${params.fdr} \
    ${input.symbol_file} \
    ${params.importfunc} \
    ${params.allele_info} \
    ${params.tx2gene_file} \
    ${params.rmdTemplate}

Replace each variable by the corresponding required file. The required files are indicated in the DESeq2 rule.

Updating/adding new tools to the workflows
-----------------------------------------------

Several yaml files provided under the folder ``shared/rules/envs`` are used to define the tools which are executed via each workflow. Here is an example from the HiC conda env::

    name: hic_conda_env_1.0
    channels:
     - conda-forge
     - anaconda
     - bioconda
    dependencies:
     - hicexplorer = 2.1.4
     - bwa = 0.7.17
     - samtools = 1.8
     - python-dateutil = 2.7.3

This file can be pointed out to the ``conda`` directive of any rule, under ``shared/rules``. Example below ::

    rule get_restrictionSite:
        input:
            genome_fasta
        output:
            enzyme + ".bed"
        params:
            res_seq = get_restriction_seq(enzyme)
        conda: CONDA_HIC_ENV
        shell:
            "findRestSite -f {input} --searchPattern {params.res_seq} -o {output} > {log.out} 2> {log.err}"

Where CONDA_HIC_ENV points to the location of the above yaml file. Under snakePipes all such global variables are defined under ``common_functions.py``

Therefore in order to change or upgrade a tool version, all you need to do is to edit the ``dependencies`` key in the yaml file to point to the new/modified tool version!


Modifying or adding new rules to the workflows
------------------------------------------------

Modifying or adding new rules to snakePipes workflows is relatively easy. Considering you want to add a new Rscript that performs a downstream analysis on the DESeq2 output in mRNAseq workflow. These would be the steps needed:

    * Test the Rscript on command line first, then move it in the ``shared/rscripts`` folder.

    * Add a rule that called the Rscript and put it under ``shared/rules`` folder.

    * Add the corresponding ``rule all``, that defines the expected output into ``workflows/mRNAseq/Snakefile``

    * Now, for easy and reproducible execution of the rule, add a ``conda`` directive and point it to the relevant conda env under ``shared/rules/envs``. Since your rule might need a new R package, `search whether it's available <https://anaconda.org/search?q=knitr>`__ in one of the conda channels and add the package name (as indicated in the conda channel) and version under the ``dependencies`` key.

    * Finally, modify the command line wrapper (``workflows/mRNAseq/mRNAseq``) to make this new feature available to the users!


Using AWS or other cloud platforms
----------------------------------

There is nothing particularly special about performing computations on AWS or other cloud platforms. Below are a few recommendations, using AWS as an example:

 1. Use a large compute node for initial installation. On AWS a ``t2.large`` node is sufficient for general installation since conda will need a couple of GB RAM for dependency resolution during setup.
 2. If you need to create custom indices, then you will need a node with at least 80GB RAM and 10 cores.
 3. Ensure that you install snakePipes on a separate EBS (or equivalent) storage block. We found that a 200GB ``/data`` partition was most convenient. This absolutely must not be the ``/`` partition, as mounting such a persistent image on other instances will result in paths being changed, which result in needing to modify large numbers of files.
 4. It's usually sufficient to use a single large (e.g., ``m5.24xlarge``) compute node, with 100+ cores and a few hundred GB RAM. This allows one to use the ``--local`` option and not have to deal with the hassle of setting up a proper cluster on AWS. Make sure the then set ``-j`` to the number of available cores on the node, so snakePipes can make the most efficient use of the resources (and minimize your bill).

Below is an example of running the mRNAseq pipeline on AWS using the resources outlined above. Note that it's best to store your input/output data on a separate storage block, since its lifetime is likely to be shorter than that of the indices.

.. code:: bash

    # Using a t2.small
    sudo mkdir /data
    mount /dev/sdf1 /data # /dev/sdf1 is a persistent storage block!
    sudo chown ec2-user /data
    cd /data

    # get datasets
    mkdir indices
    wget https://zenodo.org/record/1475957/files/GRCm38_gencode_snakePipes.tgz?download=1
    mv GRC* indices/GRCm38.tgz
    cd indices
    tar xf GRCm38.tgz
    rm GRCm38.tgz
    cd ..
    mkdir data
    wget some_data_url
    mv snakePipes_files.tar data/
    cd data
    tar xf snakePipes_files.tar
    rm snakePipes_files.tar
    cd ..

    # Edit the yaml file under indices to point to /data/indices

    # Get conda
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p conda
    export PATH=/data/conda/bin:$PATH
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    conda create -n snakePipes -c mpi-ie -c conda-forge -c bioconda snakePipes
    conda activate snakePipes
    rm Miniconda3-latest-Linux-x86_64.sh

    # setup snakePipes
    snakePipes createEnvs --only CONDA_SHARED_ENV CONDA_RNASEQ_ENV

    # Update defaults.yaml to use /data/tmp for temporary space

Then a larger instance can be spun up and the `mRNAseq` pipeline run as normal.

.. code:: bash

    mkdir /data
    mount /dev/sdf1 /data
    chown ec2-user /data
    export PATH=/data/snakePipes/bin:$PATH
    conda activate snakePipes
    mRNAseq -m alignment -i /data/data -o /data/output --local -j 192 /data/indices/GRCm28.yaml

Receiving emails upon pipeline completion
-----------------------------------------

SnakePipes can send an email to the user once a pipeline is complete if users specify ``--emailAddress``. In order for this to work, the following values need to be set in ``defaults.yaml``:

:smtpServer: The address of the outgoing SMTP server
:smtpPort: The port on the SMTP server to use (0 means to use the standard port)
:onlySSL: Set this to "True" if your SMTP server requires a full SSL connection from the beginning.
:emailSender: The name of the "user" that sends emails (e.g., snakepipes@your-domain.com)

There are two additional parameters that can be set: ``smtpUsername`` and ``smtpPassword``. These are relevant to SMTP servers that require authentication to send emails. On shared systems, it's important to ensure that other users cannot read your defaults.yaml file if it includes your password!
