__description__ = """
MPI-IE workflow for creating HiC matrices with pairtools

usage example:
    makePairs -i input-dir -o output-dir --VCFfile vcf --strains s1,s2 dm6
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
                         "downsample": False,
                         "trim": False, "trimmer": "fastp",
                         "trimmerOptions": None, "fastqc": False,
                         "reads": ["_R1", "_R2"], "ext": ".fastq.gz",
                         "fromBAM": False, "bamExt": ".bam",
                         "aligner": "bwa",
                         "alignerOptions": "-SPu -T0",
                         "plotFormat": "png",
                         "UMIDedup": False,
                         "UMIDedupOpts": "", "bcPattern": "NNNNCCCCCCCCC",
                         "UMIDedupSep": "_", "UMIBarcode": False}):
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

    parserCommon.commonOptions(optional, defaults, bw=False)

    optional.add_argument("--aligner",
                          help="Program used for mapping: BWA \
        (default: '%(default)s'). If you change this, please change \
        --alignerOptions to match.",
                          choices=["bwa"],
                          default=defaults["aligner"])

    optional.add_argument("--alignerOptions",
                          help="aligner option string, \
                          e.g.: '-SPu -T0' (default: '%(default)s')",
                          default=defaults["alignerOptions"])

    optional.add_argument("--fromBAM",
                          action="store_true",
                          help="Input folder with bam files. If provided, \
        the analysis will start from this point. If bam files contain single \
        ends, please specify --singleEnd additionally.",
                          default=defaults["fromBAM"])

    return parser


def main():
    baseDir, workflowDir, defaults = cf.setDefaults(os.path.basename(__file__))
    parser = parse_args(defaults)
    args = parser.parse_args()
    args, defaults = cf.handleUserArgs(args, defaults, parse_args)

    # add baseDir to config, although we don't use them in the Snakefile
    args.baseDir = baseDir

    # Common arguments
    cf.checkCommonArguments(args, baseDir, outDir=True)
    args.VCFfile = os.path.abspath(args.VCFfile)

    if args.fromBAM:
        args.aligner = "EXTERNAL_BAM"

#     ## End workflow-specific checks

    # Handle YAML and log files
    snakemake_cmd = \
        cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    # CreateDAG (with --DAG flag)
    cf.plot_DAG(args, snakemake_cmd, __file__, defaults)


if __name__ == "__main__":
    main()
