#!/usr/bin/env python3

__description__ = """
MPI-IE workflow for Hi-C analysis

usage example:
    HiC -i input-dir -o output-dir mm10
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
                         "downsample": False, "trim": False,
                         "trimmer": "cutadapt", "trimmerOptions": "",
                         "fastqc": False, "aligner": None, "binSize": 10000,
                         "noTAD": False,
                         "RFResolution": False, "correctionMethod": "KR",
                         "enzyme": "HindIII", "restrictRegion": None,
                         "mergeSamples": False, "nBinsToMerge": 0,
                         "findTADParams": '--thresholdComparisons 0.01',
                         "reads": ["_R1", "_R2"], "ext": ".fastq.gz",
                         "noCorrect": False, "distVsCount": False,
                         "distVsCountParams": None, "UMIDedup": False,
                         "UMIDedupOpts": "", "bcPattern": "NNNNCCCCCCCCC",
                         "UMIDedupSep": "_", "UMIBarcode": False}):
    """
    Parse arguments from the command line.
    """
    mainArgs = parserCommon.mainArguments(defaults)

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        parents=[mainArgs],
        add_help=False
    )

    # Workflow options
    optional = parser.add_argument_group('Options')

    parserCommon.commonOptions(optional, defaults, bw=False, plots=False)

    optional.add_argument("--aligner",
                          help="Program used for mapping: bwa or bwa-mem2 \
                               (default: '%(default)s').",
                          choices=["bwa", "bwa-mem2"],
                          default=defaults["aligner"])

    optional.add_argument("--RFResolution",
                          action="store_true",
                          help="Create Hi-C matrices at the restriction "
                               "fragment resolution. Using `RFResolution` "
                               "would override the --binSize argument. "
                               "(default: '%(default)s')",
                          default=defaults["RFResolution"])

    optional.add_argument("--enzyme",
                          choices=['DpnII', 'HindIII'],
                          help="Which enzyme was used to create Hi-C "
                               "library (default: '%(default)s')",
                          default=defaults["enzyme"])

    optional.add_argument("--binSize",
                          type=int,
                          metavar="INT",
                          help="Create Hi-C matrices at the given binSize. "
                               "This option is mutally exclusive with the "
                               "`--RFResolution` option \
                               (default: '%(default)s')",
                          default=defaults["binSize"])

    optional.add_argument("--restrictRegion",
                          type=str,
                          metavar="STR",
                          help="Restrict building of HiC Matrix to given "
                               "region [Chr:Start-End]. "
                               "Only one chromosome can also be specified "
                               "(default: '%(default)s')",
                          default=defaults["restrictRegion"])

    optional.add_argument("--mergeSamples",
                          action="store_true",
                          help="Merge HiC matrices and create a new matrix."
                               " If this option is specified togather with "
                               "`--sampleInfo` (see below), the samples would "
                               "be merged based on the defined groups. "
                               "(default: '%(default)s')",
                          default=defaults["mergeSamples"])

    optional.add_argument("--nBinsToMerge",
                          type=int,
                          help="If > 0 , create a lower resolution HiC matrix "
                               "for each sample by merging the given number "
                               "of bins. (default: '%(default)s')",
                          default=defaults["nBinsToMerge"])

    optional.add_argument("--findTADParams",
                          type=str,
                          metavar="STR",
                          help="parameters for HiCFindTADs. \
                               (default: '%(default)s')",
                          default=defaults["findTADParams"])

    optional.add_argument("--noTAD",
                          action="store_true",
                          help="Stop the pipeline before TAD calling. \
                               (default: '%(default)s')",
                          default=defaults["noTAD"])

    optional.add_argument("--noCorrect",
                          action="store_true",
                          help="Stop the pipeline before ICE-correction "
                               "(i.e. run only upto building the matrix). "
                               "(default: '%(default)s')",
                          default=defaults["noCorrect"])

    optional.add_argument("--distVsCount",
                          action="store_true",
                          help="Produce a plot of the ICE-corrected HiC "
                               "counts as a function of distance. This plot "
                               "could be used for QC as well as comparison "
                               "between samples for biological effects. "
                               "The plot is create using the tool "
                               "'hicDistVsCount'. (default: '%(default)s')",
                          default=defaults["distVsCount"])

    optional.add_argument("--distVsCountParams",
                          type=str,
                          metavar="STR",
                          help="parameters to run hicDistVsCount. "
                          "(default: '%(default)s')",
                          default=defaults["distVsCountParams"])

    optional.add_argument("--sampleSheet",
                          dest="sampleSheet",
                          type=str,
                          metavar="STR",
                          help="A .tsv file containing sample names and "
                               "associated groups. If provided, the file "
                               "would be used to identify groups to merge "
                               "the samples. An example can be found at "
                               "'docs/content/sampleSheet.example.tsv'"
                               "(default: None)")

    optional.add_argument("--correctionMethod",
                          dest="correctionMethod",
                          type=str,
                          metavar="STR",
                          help="Method to be used to balance the hic matrix."
                               " Available options are KR and ICE. "
                               "(default: '%(default)s')",
                               default=defaults["correctionMethod"])

    return parser


def main():
    baseDir, workflowDir, defaults = cf.setDefaults(os.path.basename(__file__))

    # get command line arguments
    parser = parse_args(defaults)
    args = parser.parse_args()
    args, defaults = cf.handleUserArgs(args, defaults, parse_args)

    # add these paths to config, although we don't use them in the Snakefile
    args.baseDir = baseDir

    # Common arguments
    cf.checkCommonArguments(args, baseDir, outDir=True)

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(
         baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    # CreateDAG
    cf.plot_DAG(args, snakemake_cmd, __file__, defaults)


if __name__ == "__main__":
    main()
