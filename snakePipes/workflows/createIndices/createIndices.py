__description__ = """
Create indices for use by snakePipes. A YAML file will be created by default in the default location where snakePipes looks for organism YAML files.

usage example:
    createIndices -o output-dir --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --blacklist blacklist.bed --ignoreForNormalization ignore.txt GRCm38_release93
"""


import argparse
import os
import sys
import textwrap
import snakePipes.common_functions as cf
import snakePipes.parserCommon as parserCommon


def parse_args(defaults={"configFile": None, "clusterConfigFile": None,
                         "maxJobs": 5, "snakemakeOptions": "--use-conda",
                         "tempDir": None, "verbose": False, "spikeinExt": None, "salmonIndexOptions": None, "eisaR_flank_length": None }):
    """
    Parse arguments from the command line.
    """
    mainArgs = parserCommon.mainArguments(defaults, createIndices=True)

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        parents=[mainArgs],
        add_help=False
    )

    parser.add_argument("genome", metavar="GENOME", help="The name to save this genome as. No spaces or special characters! Specifying an organism that already exists will cause the old information to be overwritten. See also the --userYAML option.")

    # Required arguments, which already exists as an argument group
    required = [grp for grp in parser._action_groups if grp.title == 'Required Arguments'][0]
    required.add_argument("--genomeURL",
                          required=True,
                          help="URL or local path to where the genome fasta file is located. The file may optionally be gzipped.")

    required.add_argument("--gtfURL",
                          help="URL or local path to where the genome annotation in GTF format is located. GFF is NOT supported. The file may optionally be gzipped. If this file is not specified, then RNA-seq related tools will NOT be usable.")

    # Workflow options
    optional = parser.add_argument_group('Options')

    optional.add_argument("--spikeinGenomeURL",
                          help="URL or local path to where the spikein genome fasta file is located. The file may optionally be gzipped.")

    optional.add_argument("--spikeinGtfURL",
                          help="URL or local path to where the spikein genome annotation in GTF format is located. GFF is NOT supported. The file may optionally be gzipped.")

    optional.add_argument("--spikeinExt",
                          dest="spikeinExt",
                          help="Extention of spikein chromosome names in the hybrid genome. (default: '%(default)s') .",
                          default=defaults["spikeinExt"])

    optional.add_argument("--tools",
                          help="Only produce indices for the following tools (by default, all indices will be created). The default is 'all'. 'none' will create everything except aligner indices.",
                          default="all",
                          nargs="+",
                          choices=['all', 'bowtie2', 'hisat2', 'bwa', 'bwa-mem2', 'bwameth', 'bwameth2', 'salmon', 'star', 'none'])

    optional.add_argument("--effectiveGenomeSize",
                          type=int,
                          help="The effective genome size. If you don't specify a value then the number of non-N bases will be used.")

    optional.add_argument("--spikeinBlacklist",
                          help="An optional URL or local path to a file to use to blacklist spikein organism regions (such as that provided by the ENCODE consortium).")

    optional.add_argument("--blacklist",
                          help="An optional URL or local path to a file to use to blacklist regions (such as that provided by the ENCODE consortium).")

    optional.add_argument("--ignoreForNormalization",
                          help="An optional file list, with one entry per line, the chromosomes to ignore during normalization. These are typically sex chromosomes, mitochondrial DNA, and unplaced contigs.")

    optional.add_argument("--rmskURL",
                          help="URL or local path to where the repeat masker output file is located. This is only required if you plan to run the ncRNAseq workflow.")

    optional.add_argument("--userYAML",
                          action="store_true",
                          help="By default, this workflow creates an organism YAML file where snakePipes will look for it by default. If this isn't desired (e.g., you don't want the organism to be selectable by default or you don't have write permissions to the snakePipes installation) you can specify this option and the YAML file will instead be created in the location specified by the `-o` option.")

    optional.add_argument("--salmonIndexOptions",
                          help="Options to pass to salmon for index creation.",
                          default=defaults["salmonIndexOptions"])

    optional.add_argument("--eisaR_flank_length",
                          help="Length by which to extend intronic regions with eisaR.",
                          default=defaults["eisaR_flank_length"])

    return parser


def main():
    baseDir, workflowDir, defaults = cf.setDefaults(os.path.basename(__file__))

    # get command line arguments
    parser = parse_args(defaults)
    args = parser.parse_args()
    args, defaults = cf.handleUserArgs(args, defaults, parse_args)

    # we also add these paths to config, although we don't use them in the Snakefile
    args.baseDir = baseDir

    # Common arguments
    cf.checkCommonArguments(args, baseDir, outDir=True, createIndices=True)


    ### Workflow-specific arguments
    if args.ignoreForNormalization:
        args.ignoreForNormalization = os.path.abspath(args.ignoreForNormalization)
        if not os.path.exists(args.ignoreForNormalization):
            sys.exit("The file specified by `--ignoreForNormalization` does not exist!\n")
    if args.blacklist:
        if os.path.exists(args.blacklist):
            args.blacklist = os.path.abspath(args.blacklist)
    ###

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

        # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.print_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
