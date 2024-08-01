#!/usr/bin/env python3

__description__ = """
MPI-IE workflow for scRNA-seq (CEL-Seq2 and related protocols)

usage example:
    scRNAseq -i input-dir -o output-dir mm10
"""


import argparse
import os
import sys
import textwrap
import snakePipes.common_functions as cf
import snakePipes.parserCommon as parserCommon


def parse_args(defaults={"verbose": False, "configFile": None,
                         "clusterConfigFile": None, "maxJobs": 5,
                         "snakemakeOptions": "--use-conda", "tempdir": None,
                         "downsample": False, "trim": False,
                         "reads": ["_R1", "_R2"], "ext": ".fastq.gz", "ext": ".fastq.gz",
                         "alignerOptions": "--outBAMsortingBinsN 30 --twopassMode Basic",
                         "filterGTF": "-v -P 'decay|pseudogene' ",
                         "bwBinSize": 10,
                         "plotFormat": "png",
			 "mode": "STARsolo", "BCwhiteList": None, "STARsoloCoords": ["1","7","8","7"], "myKit": "CellSeq384", "skipVelocyto": False,
                         "alevinLibraryType":"ISR", "prepProtocol":None, "expectCells" : None }):
    """
    Parse arguments from the command line.
    """
    mainArgs = parserCommon.mainArguments(defaults, workingDir=False)

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        parents=[mainArgs],
        add_help=False
    )

    # Workflow options
    optional = parser.add_argument_group('Options')

    optional.add_argument("--mode",
                          dest="mode",
                          metavar="STR",
                          help="Analysis mode. Possible settings are 'Gruen, STARsolo and Alevin' Default: %(default)s",
			  choices=['STARsolo','Alevin'],
                          type=str,
                          default=defaults["mode"])

    optional.add_argument("--downsample",
                          dest="downsample",
                          metavar="INT",
                          help="Downsample the given number of reads randomly from of each FASTQ file",
                          type=int,
                          default=defaults["downsample"])


    optional.add_argument("--filterGTF",
                          dest="filterGTF",
                          type=str,
                          help="filter annotation GTF by grep for feature counting, e.g. use --filterGTF='-v pseudogene'; "
                          "(default: '%(default)s')",
                          default=defaults["filterGTF"])


    optional.add_argument("--BCwhiteList",
                          metavar="STR",
                          help="Path to a one-column txt file with barcode whitelist. Required for the STARsolo mode,optional for Alevin mode. (default: '%(default)s')",
                          default=defaults["BCwhiteList"])

    optional.add_argument("--STARsoloCoords",
                          type=str,
                          help="Comma-separated list of values: UMI start position, UMI length, CB start position, CB length. Required for the STARsolo mode (default: '%(default)s')",
                          default=defaults["STARsoloCoords"])

    optional.add_argument("--bwBinSize",
                          dest="bwBinSize",
                          help="Bin size of output files in bigWig format (default: '%(default)s')",
                          type=int,
                          default=defaults["bwBinSize"])

    optional.add_argument("--plotFormat",
                          dest="plotFormat",
                          choices=['png', 'pdf', 'None'],
                          metavar="STR",
                          type=str,
                          help="Format of the output plots from deeptools. Select 'none' for no plot (default: '%(default)s')",
                          default=defaults["plotFormat"]),


    optional.add_argument("--myKit",
                          choices=['10Xv2', '10Xv3', 'CellSeq192', 'CellSeq384', 'Custom'],
                          metavar="STR",
                          type=str,
                          help="Library preparation kit and version to use preset barcode whitelist and CB/UMI positions for (default: '%(default)s')",
                          default=defaults["myKit"])

    optional.add_argument("--skipVelocyto",
                          dest="skipVelocyto",
                          action="store_true",
                          default=defaults["skipVelocyto"],
                          help="Skip bam filtering and generating RNA velocity counts by velocyto to save time and memory usage. (default: '%(default)s')")

    optional.add_argument("--prepProtocol",
                          dest="prepProtocol",
                          choices=["dropseq","chromiumV3", "chromium", "gemcode","citeseq","celseq","celseq2","quartzseq2"],
                          default=defaults["prepProtocol"],
                          help="Alevin mode. Specify the library prep method. (default: '%(default)s')")


    optional.add_argument("--alevinLibraryType",
                          dest="alevinLibraryType",
                          choices=["ISR", "ISF", "MSF", "MSR", "OSR", "OSF"],
                          default=defaults["alevinLibraryType"],
                          help="Alevin mode. Library orientation type. (default: '%(default)s')")

    optional.add_argument("--expectCells",
                          type=int,
                          dest="expectCells",
                          default=defaults["expectCells"],
                          help="Alevin mode. Optional to fill in if you know how many cells are expected. (default: '%(default)s')")


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

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.print_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
