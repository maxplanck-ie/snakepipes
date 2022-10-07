#!/usr/bin/env python3

__description__ = """
MPI-IE workflow for WGBS analysis

usage example:
    WGBS -i read_input_dir -o output-dir mm10
"""

import argparse
import os
import sys
import textwrap
import snakePipes.common_functions as cf
import snakePipes.parserCommon as parserCommon


def minOne(x):
    """
    Ensure x is an integer >=1
    """
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("Minimum value of 1 for --maxDist and --minCpGs")
    return x


def bounded01(x):
    """
    Ensure values are bounded [0, 1)
    """
    x = float(x)
    if x < 0 or x >= 1:
        raise argparse.ArgumentTypeError("The values for --minMethDiff and --FDR must be between 0 and 1")
    return x


def parse_args(defaults={"verbose": False, "configFile": None,
                         "clusterConfigFile":None, "maxJobs": 12,
                         "snakemakeOptions": "--use-conda", "tempDir": None,
                         "downsample": None, "aligner": None, "trim": True, "trimmer": "fastp",
                         "trimmerOptions": "-q 5 -l 30 -M 5", "maxDist": 300,
                         "minCpGs": 10, "minMethDiff": 0.1, "FDR": 0.1,
                         "MethylDackelOptions": "--mergeContext --maxVariantFrac 0.25 --minDepth 4",
                         "fastqc": False, 'bwBinSize': 25, 'plotFormat': 'png',
                         "bamExt": '.bam', "reads": ["_R1", "_R2"], "ext": ".fastq.gz",
                         "DMRprograms": "metilene,dmrseq", "minCoverage": 5,
                         "blacklist": None, "sampleSheet": None,
                         "noAutoMethylationBias": False, "fromBAM": False, "skipDOC": False,
                         "UMIDedup": False,
                         "UMIDedupOpts": "", "bcPattern": "",
                         "UMIDedupSep": "_", "UMIBarcode": False, "targetRegions": None, "metileneOptions": None, "skipBamQC": False }):

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

    parserCommon.commonOptions(optional, defaults)

    optional.add_argument("--aligner",
                          help="Program used for mapping: bwameth or bwameth2 (default: '%(default)s'). Specifying bwameth2 will run bwameth with bwa-mem2 underneath.",
                          choices=["bwameth","bwameth2"],
                          default=defaults["aligner"])

    optional.add_argument("--blacklist",
                          dest="blacklist",
                          help="Bed file(s) with positions to mask for methylation calling. Useful for masking SNPs in your strain of interest. (default: '%(default)s')",
                          default=defaults["blacklist"])

    optional.add_argument("--targetRegions",
                          dest="targetRegions",
                          help="Bed file(s) with regions of interest to evaluate methylation for. (default: '%(default)s')",
                          default=defaults["targetRegions"])

    optional.add_argument("--sampleSheet",
                          dest="sampleSheet",
                          help="Perform differential methylation analysis between groups of samples by providing a text file with sample information to use for statistical analysis. (default: '%(default)s')",
                          default=defaults["sampleSheet"])

    optional.add_argument("--noAutoMethylationBias",
                          action="store_true",
                          default=defaults["noAutoMethylationBias"],
                          help="If specified, MethylDackel mbias will NOT be run and the suggested parameters from it will NOT be used for methylation extraction. You can instead supply them manually in --MethylDackelOptions.")

    optional.add_argument("--maxDist",
                          type=minOne,
                          help="The maximum distance between CpGs in a DMR (for metilene, default: '%(default)s')",
                          default=defaults["maxDist"])

    optional.add_argument("--minCpGs",
                          type=minOne,
                          help="The minimum number of CpGs in a DMR (for metilene, default: '%(default)s')",
                          default=defaults["minCpGs"])

    optional.add_argument("--minMethDiff",
                          type=float,
                          help="The minimum methylation change in methylation for CpG inclusion in DMR detection (for metilene, default: '%(default)s')",
                          default=defaults["minMethDiff"])

    optional.add_argument("--minCoverage",
                          type=int,
                          help="The minimum coverage needed for across all samples for a CpG to be used in DMR calling and PCA. Note that you can change this value without overwritting the DMR output.  Default: '%(default)i'",
                          default=defaults["minCoverage"])

    optional.add_argument("--FDR",
                          type=float,
                          help="FDR threshold for returned DMRs (default: '%(default)s')",
                          default=defaults["FDR"])

    optional.add_argument("--MethylDackelOptions",
                          help="Options to pass to MethylDackel extract. You are highly advised NOT to set a minimum coverage at this step. Default: '%(default)s'",
                          default=defaults["MethylDackelOptions"])

    optional.add_argument("--fromBAM",
                         help="If specified, the input is taking from BAM files containing alignments rather than fastq files. See also --bamExt.",
                         action="store_true")

    optional.add_argument("--skipBamQC",
                         help="If specified, filtering of input bam files, as well as calculation of quality metrics, will be skipped. ",
                         action="store_true")

    optional.add_argument("--bamExt",
                          help="If --fromBAM is specified, this is the expected file extension. Removing it yields sample names. Default: '%(default)s'",
                          default=defaults["bamExt"])

    optional.add_argument("--singleEnd",
                          dest="pairedEnd",
                          help="If --fromBAM is specified, this indicates that the input BAM files contain paired-end data. The option is ignored unless --fromBAM is given.",
                          action="store_false")

    optional.add_argument("--DMRprograms",
                          help="If a sample sheet is provided, use the specified DMR-calling programs. Multiple programs can be comma-separated with no spaces (e.g., 'metilene,dmrseq,DSS'). The available programs are metilene, dmrseq, and DSS (note that this is very slow). Default: %(default)s.",
                          default=defaults["DMRprograms"])

    optional.add_argument("--metileneOptions",
                          help="Options to pass to metilene. Default: '%(default)s'",
                          default=defaults["metileneOptions"])

    optional.add_argument("--skipDOC",
                         action="store_true",
                         help="Skip depth of coverage calculation with deepTools.",default=defaults["skipDOC"])

    optional.add_argument("--GCbias",
                         action="store_true",
                         help="Perform GC bias calculation with deepTools.")

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
