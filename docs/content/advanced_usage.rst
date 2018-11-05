Advanced usage of snakePipes
============================

Calling snakemake without using the wrapper script
--------------------------------------------------

Please save a copy of this config yaml file and provide an adjusted config via '--configfile' parameter!

example call::

    snakemake --snakefile /path/to/snakemake_workflows/workflows/ATAC-seq/Snakefile
              --configfile /path/to/snakemake_workflows/workflows/ATAC-seq/defaults.yaml
              --directory /path/to/outputdir
              --cores 32

Using AWS or other cloud platforms
----------------------------------

There is nothing particularly special about performing computations on AWS or other cloud platforms. Below are a few recommendations, using AWS as an example:

 1. Use a small compute node for initial installation. On AWS a `t2.small` node is sufficient for general installation since conda will need 1-2GB RAM for dependency resolution during setup.
 2. If you can need to create custom indices, then you will need a node with at least 80GB RAM and 10 cores.
 3. Ensure that you install snakePipes on a separate EBS (or equivalent) storage block. We found that a 200GB `/data` partition was most convenient. This absolutely must not be the `/` partition, as mounting such a persistent image on other instances will result in paths being changed, which result in needing to modify large numbers of files.
 4. It's usually sufficient to use a single large (e.g., `m5.24xlarge`) compute node, with 100+ cores and a few hundred GB RAM. This allows one to use the `--local` option and not have to deal with the hassle of setting up a proper cluster on AWS. Make sure the then set `-j` to the number of available cores on the node, so snakePipes can make the most efficient use of the resources (and minimize your bill).

Below is an example of running the RNA-seq pipeline on AWS using the resources outlined above. Note that it's best to store your input/output data on a separate storage block, since its lifetime is likely to be shorter than that of the indices.

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
    source activate snakePipes
    rm Miniconda3-latest-Linux-x86_64.sh

    # setup snakePipes
    snakePipes createEnvs --only CONDA_SHARED_ENV CONDA_RNASEQ_ENV

    # Update defaults.yaml to use /data/tmp for temporary space

Then a larger instance can be spun up and the `RNA-seq` pipeline run as normal.

.. code:: bash

    mkdir /data
    mount /dev/sdf1 /data
    chown ec2-user /data
    export PATH=/data/snakePipes/bin:$PATH
    source activate snakePipes
    RNA-seq -m alignment -i /data/data -o /data/output --local -j 192 /data/indices/GRCm28.yaml
