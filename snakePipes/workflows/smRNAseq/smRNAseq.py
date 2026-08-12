__description__ = """
MPI-IE workflow for RNA mapping and analysis

usage example:
    smRNAseq -i input-dir -o output-dir mm10
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
                         "mode": "alignment,deepTools_qc", "downsample": False,
                         "trim": True, "trimmer": "fastp",
                         "fastqc": False,
                         "libraryType": 2,
                         "sampleSheet": None,
                         "formula": "",
                         "reads": ["_R1", "_R2"], "ext": ".fastq.gz",
                         "bwBinSize": 25, "plotFormat": "png",
                         "pairedEnd": True,
                         "fdr": 0.05,
                         "tesmallMinLen": 16, "tesmallMaxLen": 36,
                         "tesmallOptions": None,
                         "fromBAM": False}):

    """
    Parse arguments from the command line.
    """
    mainArgs = parserCommon.mainArguments(defaults, workingDir=False)
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
                          help="workflow running modes (available: alignment, deepTools_qc)"
                          " (default: '%(default)s')",
                          default=defaults["mode"])

#    parserCommon.commonOptions(optional, defaults, bw=True)

    optional.add_argument("--libraryType",
                          help="user provided library type strand specificity. featureCounts style: 0, 1, 2 (Illumina TruSeq); default: '%(default)s')",
                          type=int,
                          default=defaults["libraryType"])

    optional.add_argument("--sampleSheet",
                          help="Information on samples (required for DE analysis); see "
                               "'https://github.com/maxplanck-ie/snakepipes/tree/master/docs/content/sampleSheet.example.tsv' for example."
                               " The column names in the tsv files are 'name' and 'condition'. The first entry"
                               " defines which group of samples are control. "
                               " This way, the order of comparison and likewise the sign of values can be changed."
                               " The DE analysis might fail if your sample names begin with a number. So watch out"
                               " for that! (default: '%(default)s')",
                          default=defaults["sampleSheet"])

    optional.add_argument("--formula",
                          dest="formula",
                          help="Design formula to use in linear model fit (default: '%(default)s')",
                          default=defaults["formula"])

    optional.add_argument("--fromBAM",
                         action="store_true",
                         help="Input folder with bam files. If provided, the analysis will start from this point. If bam files contain single ends, please specify --singleEnd additionally.",
                         default=defaults["fromBAM"])

    optional.add_argument("--singleEnd",
                          dest="pairedEnd",
                          action="store_false",
                          help="input data is single-end, not paired-end. This is only used if --fromBAM is specified.")

    optional.add_argument("--FDR",
                          dest="fdr",
                          help="FDR threshold to apply for filtering DE genes"
                               "(default: '%(default)s')",
                          default=defaults["fdr"])

    optional.add_argument("--tesmallMinLen",
                          type=int,
                          help="Passed to TEsmall as -m/--minlen: discard trimmed reads shorter "
                               "than this (reads too short even before adapter removal are also "
                               "discarded). (default: %(default)s)",
                          default=defaults["tesmallMinLen"])

    optional.add_argument("--tesmallMaxLen",
                          type=int,
                          help="Passed to TEsmall as -M/--maxlen: discard trimmed reads longer "
                               "than this (reads too long even before adapter removal are also "
                               "discarded). (default: %(default)s)",
                          default=defaults["tesmallMaxLen"])

    optional.add_argument("--tesmallOptions",
                          help="Additional TEsmall option string, appended as-is to the TEsmall "
                               "call, e.g.: '--maxaln 200 --mismatch 1'. Use this for any TEsmall "
                               "argument not already covered by its own snakePipes flag (see "
                               "'TEsmall --help' for the full list). (default: '%(default)s')",
                          default=defaults["tesmallOptions"])

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

    ## Begin workflow-specific checks
    modeTemp = args.mode.split(",")
    validModes = set(["deepTools_qc"])
    for mode in modeTemp:
        if mode not in validModes:
            sys.exit("{} is not a valid mode!\n".format(mode))

    if args.tesmallMinLen >= args.tesmallMaxLen:
        sys.exit("--tesmallMinLen ({}) must be smaller than --tesmallMaxLen ({})!\n".format(
            args.tesmallMinLen, args.tesmallMaxLen))

    ## End workflow-specific checks

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.plot_DAG(args,snakemake_cmd, __file__,defaults)
