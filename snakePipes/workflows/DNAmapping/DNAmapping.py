__description__ = """
MPI-IE workflow for DNA mapping

usage example:
    DNAmapping -i input-dir -o output-dir mm10
"""


import argparse
import os
import sys
import textwrap
import snakePipes.common_functions as cf
import snakePipes.parserCommon as parserCommon


def parse_args(defaults={"verbose": False, "configFile": None,
                         "clusterConfigFile": None, "maxJobs": 5,
                         "snakemakeOptions": "--use-conda", "tempDir": None,
                         "mode": "mapping", "downsample": False, "trim": False,
                         "trimmer": "cutadapt", "trimmerOptions": "", "fastqc": False,
                         "qualimap": False, "dedup": False, "ext": ".fastq.gz",
                         "properPairs": False, "insertSizeMax": 1000,
                         "GCBias": False, "reads": ["_R1", "_R2"],
                         "bwBinSize": 25, "mapq": 0, "plotFormat": "png",
                         "alignerOpts": "", "mateOrientation": "--fr",
                         "UMIDedup": False, "UMIDedupOpts": "",
                         "UMIDedupSep": "_", "UMIBarcode": False, "cutntag": False,
                         "bcPattern": "NNNNCCCCCCCC", "aligner":"Bowtie2"}):
    """
    Parse arguments from the command line.
    """
    mainArgs = parserCommon.mainArguments(defaults)
    snpArgs = parserCommon.snpArguments(defaults)

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        parents=[mainArgs, snpArgs],
        add_help=False
    )


    # Workflow options
    optional = parser.add_argument_group('Options')
    optional.add_argument("-m", "--mode",
                          dest="mode",
                          help="workflow running modes (available: 'mapping,"
                          "allelic-mapping')(default: '%(default)s')",
                          default=defaults["mode"])

    parserCommon.commonOptions(optional, defaults)

    optional.add_argument("--alignerOpts",
                          help="Options that will be passed to Bowtie2 or bwa. You can specify things such as `--local` or "
                          "`--very-sensitive` here. The mate orientation and maximum insert size are specified "
                          "elsewhere. Read group information is set automatically. Note that you may need to escape "
                          r"the first - (e.g., '\--very-fast'). Default: '%(default)s'.",
                          default=defaults["alignerOpts"])

    optional.add_argument("--cutntag",
                          help="if set, Bowti2 is used for mapping with parameters as has been used "
                          "in the method section of Kaya-okur et al. 2019. ('--local --very-sensitive-local "
                          "--no-mixed --no-discordant --phred33 -I 10 -X 700')"
                          "Setting this flag overwrites the '--alignerOpts' and '--insertSizeMax'."
                          " Default is '%(default)s'.",
                          action="store_true")

    optional.add_argument("--mateOrientation",
                          help="The --fr, --ff, or --rf option for bowtie2 (default: '%(default)s')",
                          default=defaults["mateOrientation"])

    optional.add_argument("--qualimap",
                          dest="qualimap",
                          action="store_true",
                          help="activate Qualimap (default: '%(default)s')",
                          default=defaults["qualimap"])

    optional.add_argument("--dedup",
                          dest="dedup",
                          action="store_true",
                          help="retain only de-duplicated reads/read pairs "
                          "(given single-/paired-end data), recommended for "
                          "ChIPseq data (default: '%(default)s')",
                          default=defaults["dedup"])

    optional.add_argument("--properPairs",
                          action="store_true",
                          help="retain only reads mapping in proper pairs (default: '%(default)s')",
                          default=defaults["properPairs"])

    optional.add_argument("--mapq",
                          dest="mapq",
                          metavar="INT",
                          help="retain only reads with at least the given "
                          "mapping quality. We recommend using"
                          "mapq of 3 or more for ChIPseq to remove all true "
                          "multimapping reads. (default: '%(default)s')",
                          type=int,
                          default=defaults["mapq"])

    optional.add_argument("--insertSizeMax",
                          help="Maximum insert size allowed during mapping (default: '%(default)s')",
                          type=int,
                          default=defaults["insertSizeMax"])

    optional.add_argument("--GCBias",
                          action="store_true",
                          help="run computeGCBias quality control "
                          "(long runtime!). Note that GCBias analysis is "
                          "skipped if downsampling is specified "
                          "(default: '%(default)s')",
                          default=defaults["GCBias"])

    optional.add_argument("--aligner",
                          help="Program used for mapping: Bowtie2 or bwa (default: '%(default)s').",
                          choices=["Bowtie2","bwa","bwa-mem2"],
                          default=defaults["aligner"])

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
    cf.checkCommonArguments(args, baseDir, outDir=True)

    ## Begin workflow-specific code
    # check for Allele-specific mapping mode
    args.allele_mode = cf.checkAlleleParams(args)
    # convert file path to abspath
    if args.allele_mode == "create_and_map":
        args.VCFfile = os.path.abspath(args.VCFfile)
    elif args.allele_mode == "map_only":
        args.SNPfile = os.path.abspath(args.SNPfile)
        args.NMaskedIndex = os.path.abspath(args.NMaskedIndex)

    ## End workflow-specific clode

    # Handle YAML and log files

    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.plot_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
