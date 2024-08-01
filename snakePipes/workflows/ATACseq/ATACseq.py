__description__ = """
MPI-IE workflow for ATACseq Analysis

usage example:
    ATACseq -d working-dir mm10
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
                         "reads": ["_R1", "_R2"],
                         "maxFragmentSize": 150, "minFragmentSize": 0,
                         "fromBAM": False ,"qval": 0.001 ,"sampleSheet": None,
                         "externalBed": None,
                         "bamExt": "filtered.bam", "fdr": 0.05,
                         "absBestLFC": 1, "UMIDedup": False,
                         "UMIDedupOpts": "", "UMIDedupSep": "_",
                         "peakCaller": "MACS2",
                         "UMIBarcode": False, "bcPattern": ""}):

    """
    Parse arguments from the command line.
    """
    mainArgs = parserCommon.mainArguments(defaults, workingDir=True)

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        parents=[mainArgs],
        add_help=False
    )



    # Workflow options
    optional = parser.add_argument_group('Options')
    optional.add_argument("--peakCaller",
                          help="The peak caller to use. The default is %(default)s",
                          choices=["MACS2", "HMMRATAC", "Genrich"],
                          default=defaults['peakCaller'])

    optional.add_argument("--maxFragmentSize",
                          help="Maximum size of (typically nucleosomal) fragments for inclusion in the analysis (default: '%(default)s')",
                          type=int,
                          default=defaults['maxFragmentSize'])

    optional.add_argument("--minFragmentSize",
                          help="Minimum size of (typically nucleosomal) fragments for inclusion in the analysis (default: '%(default)s')",
                          type=int,
                          default=defaults['minFragmentSize'])

    optional.add_argument("--qval",
                          dest="qval",
                          metavar="INT",
                          help="qvalue threshold for MACS2(default: '%(default)s')",
                          type=float,
                          default=defaults['qval'])


    optional.add_argument("--sampleSheet",
                          dest="sampleSheet",
                          help="Invoke differential accessibility analysis by providing information on samples; see 'https://github.com/maxplanck-ie/snakepipes/tree/master/docs/content/sampleSheet.example.tsv' for example. IMPORTANT: The first entry defines which group of samples are control. With this, the order of comparison and likewise the sign of values can be changed! Also, the condition `control` should not be used (reserved to mark input samples in the ChIPSeq workflow (default: '%(default)s').",
                          default=defaults["sampleSheet"])

    optional.add_argument("--externalBed",
                          help="A bed file with intervals to be tested for differential binding. (default: '%(default)s')",
                          default=defaults["externalBed"])

    optional.add_argument("--fromBAM",
                          dest="fromBAM",
                          help="Input folder with bam files. If provided, the analysis will start from this point. (default: '%(default)s')",
                          default=defaults["fromBAM"])

    optional.add_argument("--bamExt",
                         dest="bamExt",
                         help="Extention of provided bam files, will be substracted from basenames to obtain sample names. (default: '%(default)s')",
                         default=defaults["bamExt"])

    optional.add_argument("--FDR",
                          dest="fdr",
                          help="FDR threshold to apply for filtering DB regions"
                               "(default: '%(default)s')",
                          default=defaults["fdr"])

    optional.add_argument("--LFC",
                          dest="absBestLFC",
                          help="Log fold change threshold to apply for filtering DB regions"
                               "(default: '%(default)s')",
                          default=defaults["absBestLFC"])

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
    cf.checkCommonArguments(args, baseDir)

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.print_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
