.. _running_snakepipes:

Running snakePipes
==================

Pipelines under snakePipes are designed in a way such that all workflows are configured and ran in a similar way.


An example with ChIPseq data
------------------------------

A **typical ChIPseq analysis** of human samples starts from paired-end FASTQ files in the directory ``input-dir``:

.. code:: bash

    $ ls /path/to/input-dir/
    my_H3K27ac_sample_R1.fastq.gz  my_H3K27me3_sample_R1.fastq.gz  my_Input_sample_R1.fastq.gz
    my_H3K27ac_sample_R2.fastq.gz  my_H3K27me3_sample_R2.fastq.gz  my_Input_sample_R2.fastq.gz

The :ref:`ChIPseq` workflow requires the files to be processed via the :ref:`DNAmapping` workflow first. We therefore run the DNAmapping workflow :

.. code:: bash

    $ DNAmapping -i /path/to/input-dir -o /path/to/output-dir --mapq 5 -j 10 --dedup hs37d5

* ``--mapq 5`` would filter mapped reads for a minimum mapping quality of 5. This would keep only primary alignments from bowtie2, sufficient for downstream analysis.

* ``--dedup`` would remove PCR duplicates (reads with matching 5' position in the genome), a typical step in ChIPseq analysis.

* ``-j 10`` defines 10 jobs to be run in parallel on the cluster (see below).

* ``hs37d5`` is the name of the genome (keyword for the yaml). The yaml file corresponding to this genome should exist as ``snakePipes/shared/organisms/hs37d5.yaml``. (see :ref:`setting_up` for details).

All individual jobs of the workflow will be submitted to the Grid engine using the command specified under ``shared/cluster.yaml``. The parameter ``-j`` defines the number of jobs to be run in parallel, while the number of threads per job is hard-coded in the workflows.

**To run the workflow locally**, use the parameter ``--local`` for local mode and the parameter ``-j 10`` to specify the maximal number of used CPU threads (here: 10).

**For single-end FASTQ files**, Note that single end data still needs a valid suffix (e.g. sample1_R1.fastq.gz). With a proper suffix, single end mode is detected by default. When executing some workflows with the ``--fromBAM`` flag, it is still necessary to set ``--singleEnd``.

Once the DNAmapping run is finished sucessfully. We can run the ChIPseq analysis in the same directory.

.. code:: bash

    $  ChIPseq -d /path/to/dnamapping-output/ hs37d5 chip-samples.yaml

* ``-d`` specifies the directory where the output of DNAmapping workflow lies. The ChIPseq workflow would also write it's output there.

* ``hs37d5`` is the name of the genome (keyword for the yaml).

* ``chip-samples.yaml`` is a yaml file that defines for each ChIP sample, the corresponding control (input) sample and the type of mark (broad/sharp). See :ref:`ChIPseq` for more details on how to setup this yaml file.

The ChIPseq workflow would follow up from the DNAmapping outputs and perform peak calling, create ChIP-input normalized coverage files and also perform differential (control-test) analysis if a sample information file is provided (see below).

.. _sampleinfo:

The sample sheet
----------------

Most of the workflows allow users to perform grouped operations as an option, for example
differential expression analysis in mRNAseq workflow, differential binding analysis in
ChIPseq workflow, differential open-chromatin analysis in ATACseq workflow or merging of
groups in Hi-C workflow. For all this analysis, snakePipes needs a ``sampleSheet.tsv`` file (file name is not important, but it has to be tab-separated) that contains sample grouping information. In most cases users would want to groups samples by replicates. The format of the file is as follows:

::

    name   condition
    sample1    control
    sample1    control
    sample2    test
    sample2    test

The name section referes to sample names (without the read suffix), while the condition
section refers to sample group (control/test, male/female, normal/diseased etc..)

Using BAM input
---------------

In many workflows it is possible to directly use BAM files as input by specifying ``--fromBAM``. Note that you must then specify whether you have paired-end (the default) or single-end data. This is typically done with the ``--singleEnd`` option.

Changing read extensions or mate designators
--------------------------------------------

The default file names produced by Illumina sequencers are of the form ``<sample>_R1.fastq.gz`` and ``<sample_R2.fastq.gz``. However, sometimes public datasets will instead have a ``.fq.gz`` suffix or use ``_1`` and ``_2`` as mate designators. To enable this, the ``--ext`` option can be used to change ``.fastq.gz`` default suffix to ``.fq.gz`` and ``--reads`` to ``_1 _2``.

Common considerations for all workflows
----------------------------------------

All of the snakePipes workflows that begin with a FASTQ file, perform the same pre-processing steps.

* **Linking/downsampling the FASTQ file** : The FASTQ rule in the workflows links the input FASTQ file into the FASTQ folder in the output directory. If ``downsampling`` is specified, the FASTQ folder would contain the downsampled FASTQ file.

.. note:: The DNAmapping and RNA-mapping pipelines can take either single, or paired-end FASTQ files. For paired-end data, the reads ``R1`` and ``R2`` are expected to have the suffix ``_R1`` and ``_R2`` respectively, which can be modified in the ``defaults.yaml`` file using the ``reads`` key, to your needs. For example, files downloaded from NCBI would normally have the extention ``.1.fastq.gz`` and ``.2.fastq.gz``. Also, please check the ``ext`` key in the configuration file if you wish to modify the read extension (default is ``.fastq.gz``).


* **Quality/adapter trimming** (optional): If ``--trim`` is selected, the ``trimming`` rule would run the selected program (either `Trimgalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`__, or `Cutadapt <https://journal.embnet.org/index.php/embnetjournal/article/view/200/479>`__) on the files in the FASTQ folder, and would produce another folder with name ``FASTQ_<program>``, where <program> is either ``Cutadapt`` or ``Trimgalore``.

* **FastQC** (optional): If ``--fastqc`` is specified, the ``FASTQC`` rule would run `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__ on the input files and store the output under ``FastQC`` folder. If trimming is specified, FastQC is always produced on trimmed files, and stored under ``FastQC_trimmed`` folder.

* **--snakemakeOptions**: All wrappers contain a ``--snakemakeOptions`` parameter, which is quite useful as it can be used to pass on any arguments directly to snakemake. One use case is to perform a *dry run*, i.e. to check which programs would be executed and which outputs would be created by the workflow, without actually running it. This can be executed via ``--snakemakeOptions="-np"``. This would also print the commands to be used during the run.

* **--DAG**: All workflows can produce a `directed acyclic graph <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`__ of themselves, using the ``--DAG`` option in the wrappers. This could be useful in reporting/presenting the results.

* **--keepTemp**: This option control temporary/intermediate files are to be kept after the workflow is finished. Normally the temporary files are removed after analysis.

* **--bwBinSize**: This option is available for most workflows, and refers to the bin size used to create the coverage files. `BigWig files <https://genome.ucsc.edu/goldenpath/help/bigWig.html>`__ are created by most workflows in order to allow downstream analysis and visualization of outputs. This argument controls the size of the bins in which the genome is divided for creating this file. The default is sufficient for most analysis.

* **Temporary directory/files**: Some tools need additonal space during runtime (eg. ``samtools sort -T [DIR] ...``). SnakePipes uses the core tool ``mktemp`` to create temporary directories in some rules. On Linux-based systems the global env variabale ``$TMPDIR`` is honored.
  On Mac OS and if $TMPDIR is empty, we fallback to `/tmp/` as the parent temporary directory. For performance reasons, it is recommended that the $TMPDIR points to a local drive (and not eg. an NFS share). Please make sure there is enough space! 

Logging of outputs
~~~~~~~~~~~~~~~~~~~

snakePipes produces logs at three diferrent levels.

* **<workflow>.log**: This file would be generated on the working directory, and contains everything printed on the screen via snakemake and python wrappers.

* **<workflow>_organism.yaml**: This file is a copy of the YAML file specifying where all of the genomic indices, annotations, and other files are located.

* **cluster_logs**: In case snakePipes is setup with a cluster, the folder ``cluster_logs`` would contain the output and error messages from the cluster scheduler.

* **<output>/logs**: Each output folder from snakePipes workflows contain their own log (``.err`` and ``.out``) file under ``/logs/`` folder. This contains the messages directly from the executed tools.

.. note:: For most cases where a tool fails, these files contain useful debugging information. However sometimes, the error can't be captured in these files and therefore ends up in the ``cluster_logs`` folder.

Quality-Checks
~~~~~~~~~~~~~~~~~~~~~

All workflows under snakePipes employ various quality-checks (QC) to inform users of the data quality.

* **MultiQC** : All workflows in snakePipes output a ``MultiQC`` folder, which summerizes the QC metrics obtained from various tools in the workflow via `MultiQC <https://multiqc.info/>`__, in an interactive HTML report. This output is quite useful to compare samples and get an overview of the data quality from all samples.

* **deepTools**: `deepTools <deeptools.readthedocs.io>`__ are a popular set of tools that perform QC, normalization and visualization of NGS data. In snakePipes, most workflows (except HiC and scRNAseq) contain outputs from various deepTools modules on the samples. The coverage files (bigWigs), are also generated by deepTools (bamCoverage and bamCompare modules). Therefore, it's useful to look at the deepTools documentation before inspecting these results.

.. note:: We strongly encourage users to understand these quality matrices and inspect the results from QC, before making biological conclusions or preceeding to downstream analysis.
