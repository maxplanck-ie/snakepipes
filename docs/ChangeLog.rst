Change Log
================

The history of development of the workflows is listed below, along with github IDs of
people with most (but not all) contribution to the changes.


version 0.0.1 - March 23, 2016 - @kilpert
- initial version

version 0.1.0 - June 15, 2016 - @asrichter
- added --fastqc and --bw-binsize parameters to DNA-mapping wrapper script
- additional organisms are now supported by adding new genome.py files
- defined (effective) genome size as (genome length)-(number of 'N's) in genome.py files
- simplified cluster submission by taking 'threads' parameter from rule definition, removed cluster.yaml
- added DNA-mapping example configuration yaml file
- added function get_fragment_length() to DNA-mapping internals.snakefile to parse median insert size from Picard output
- added deepTools_qc.snakefile as one common snakefile for all deepTools rules, moved all deepTools rules from ChIP-seq Snakefile and computeGCBias.snakefile
- moved include statements of module snakefiles from rules.snakefile to DNA-mapping Snakefile to simplify structure, removed rules.snakefile
- renamed functions do_TrimGalore(trim) and do_InsertSizeMetrics(paired) in DNA-mapping Snakefile
- renamed load_organisms.snakefile to load_organism_data.py
- removed debugging.snakefile as it was not used by any pipeline
- removed --local-cores parameter from all wrapper scripts as there are no local snakemake rules defined
- many small other changes

version 0.1.1 - June 17, 2016 -  @asrichter
- added option to run workflow locally instead of cluster submission
- FASTQ.snakefile replaces FASTQ_symlink.snakefile and FASTQ_downsample.snakefile
- several small changes

version 0.2.0 - June 22, 2016 -  @asrichter
- added filtering option to filter BAM files for duplication, proper pairs and MAPQ
- added variable 'outdir' to configuration
- many small changes

version 0.3.0 - June 24, 2016 -  @asrichter
- rewrote ChIP-seq workflow completely including wrapper script ChIP-seq
- added histoneHMM for calling broadly enriched regions
- added MACS2 peak quality controls
- revised example.config.yaml files
- moved function get_fragment_length() to common_functions.py
- positional instead of required optional command line arguments
- many small changes

version 0.3.1 - June 25, 2016 -  @asrichter
- run Picard quality control on unfiltered BAM files
- added --gcbias parameter to DNA-mapping wrapper script to run computeGCBias optionally
- replaced --input-dir and --output-dir by --working-dir parameter in ChIP-seq wrapper script to specify the working directory, which is output directory of the pipeline and must also contain the DNA-mapping pipeline output files
- bugfixes

version 0.3.2 - June 27, 2016 -  @asrichter
- added generation of QC reports for all samples to ChIP-seq pipeline
- added consistency check for ChIP-seq pipeline whether all required input files exist for all samples
- added peak count to MACS2 peak quality controls

version 0.3.2.1 - June 27, 2016 - Andreas Richter
- added documentation to README.md
- moved R library

version 0.4 - 2016 - @kilpert, @steffenheyne

version 0.5 - 2017 - @steffenheyne, @kilpert, @mirax87
- major cleanup and refactoring of wrappers and code structure (but not rules)
- scRNAseq workflow added
- using yaml config files all over, ie.

  1) all wrappers write out a config yaml
  2) Snakefile configuration only depends on provided '--configfile'

- better handling and and more usage of common_functions
- all genome config files converted to yaml, function to load genome information
- own genome file can be provided instead of only pre-configured ones
- paths config file converted to yaml
- eliminate organism specific t2g and symbol files, they are now created by rule from gtf
- check that sleuth is only run if replicates are available
- bugfix: only provided wrapper arguments overwrite either defaults or configfile options, before
  a not provided wrapper option (but in config file provided) would be overwritten by default,
  now there is true hierachy: defaults->configfile->wrapper !


version 0.6 (a.k.a Tiger RattleSnake) - Sept 2017 - @vivekbhr
- MAJOR CHANGES:

  - Allele-Specific mapping : Allele-specific DNA and RNA-mapping is now possible and both ChIP-Seq and RNA-seq pipeline can handle "allele_mapping" mode.
  - Differential binding : Differential binding can be performed using CSAW, both normal and allele-specific.
  - MultiQC : MultiQC summarizes QC results for DNA-mapping and RNA-mapping workflows.
  - R package : New R package called snakediff contains all functions for differential binding and expression. The DE functions have been transformed into clean and consistent layout, with docs.
  - Deeptools commands : All deeptools commands now moved into common deeptools_cmds module, and are shared between functions.

- MINOR CHANGES:

  - Wrappers have been formatted for easier readability.
  - chr X, Y, M and unmapped scaffolds are now ignored for normalization in bamCoverage (DNA-mapping) and bamCompare (Chip-Seq).
  - bamCompare outputs are named with suffix <control>, instead of a static suffix "Input", this allows the workflow to run with different controls
    eg. Input, H3 and keep results in same dir. However, MACS2 would still over-write the output files from previous run (needs fix).

Version 0.6.1

- MINOR CHANGES:
  - The allele-specific option is no longer on by default (it was @vivekbhr's fault)

Version 0.6.2

- MINOR CHANGES:
  - Explicitly define which snakemake version to use

Version 0.7 (a.k.a Green Mamba) - Nov 2017 - @vivekbhr

- MAJOR CHANGES:
  - Read the Docs integration
  - New workflow Hi-C, from mapping to TAD calling, using BWA and HiCExplorer

Version 0.7.1 - Jan. 12 2018

- MINOR CHANGES:
  - Fixed issue #57 again via pull request #79, which was related to single-end ChIPseq processing
