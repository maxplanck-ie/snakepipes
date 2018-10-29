Advanced usage of snakePipes
============================

Calling snakemake without using the wrapper script:
---------------------------------------------------

It's also possible to directly run ``snakemake``. Please save a copy of this config yaml file and provide an adjusted config via ``--configfile`` parameter!

example call::

    snakemake --snakefile /path/to/snakemake_workflows/workflows/ATAC-seq/Snakefile
              --configfile /path/to/snakemake_workflows/workflows/ATAC-seq/defaults.yaml
              --directory /path/to/outputdir
              --cores 32
