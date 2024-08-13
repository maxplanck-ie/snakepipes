__description__ = """
MPI-IE workflow for ncRNAseq mapping and analysis

usage example:
    ncRNAseq -i input-dir -o output-dir mm10
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
                         "trim": False, "trimmer": "fastp",
                         "trimmerOptions": None, "fastqc": False,
                         "libraryType": 2, "aligner": "STAR",
                         "alignerOptions": "--sjdbOverhang 100 --outSAMstrandField intronMotif --outFilterMultimapNmax 1000 --outFilterMismatchNoverLmax 0.1 --outSAMattributes Standard --outSAMunmapped Within --outSAMtype BAM Unsorted",
                         "sampleSheet": None,
                         "reads": ["_R1", "_R2"], "ext": ".fastq.gz",
                         "bwBinSize": 25, "plotFormat": "png",
                         "fromBAM": False, "bamExt": ".bam", "pairedEnd": True,
                         "UMIDedup": False,
                         "UMIDedupOpts": "", "bcPattern": "NNNNCCCCCCCCC",
                         "UMIDedupSep": "_", "UMIBarcode": False}):
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
    optional.add_argument("-m", "--mode",
                          help="workflow running modes (available: 'alignment, deepTools_qc')"
                          " (default: '%(default)s')",
                          default=defaults["mode"])

    parserCommon.commonOptions(optional, defaults, bw=True)

    optional.add_argument("--aligner",
                          help="Program used for mapping: STAR (default: '%(default)s'). If you change this, please change --alignerOptions to match.",
                          default=defaults["aligner"])

    optional.add_argument("--alignerOptions",
                          help="STAR option string, e.g.: '--twopassMode Basic' (default: '%(default)s')",
                          default=defaults["alignerOptions"])

    optional.add_argument("--sampleSheet",
                          help="Information on samples (required for DE analysis); see "
                               "'https://github.com/maxplanck-ie/snakepipes/tree/master/docs/content/sampleSheet.example.tsv' for example."
                               " The column names in the tsv files are 'name' and 'condition'. The first entry"
                               " defines which group of samples are control. "
                               " This way, the order of comparison and likewise the sign of values can be changed."
                               " The DE analysis might fail if your sample names begin with a number. So watch out"
                               " for that! (default: '%(default)s')",
                          default=defaults["sampleSheet"])

    optional.add_argument("--fromBAM",
                         action="store_true",
                         help="Input folder with bam files. If provided, the analysis will start from this point. If bam files contain single ends, please specify --singleEnd additionally.",
                         default=defaults["fromBAM"])

    optional.add_argument("--bamExt",
                          help="Extention of provided bam files, will be substracted from basenames to obtain sample names. (default: '%(default)s')",
                          default=defaults["bamExt"])


    optional.add_argument("--singleEnd",
                          dest="pairedEnd",
                          action="store_false",
                          help="input data is single-end, not paired-end. This is only used if --fromBAM is specified.")

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
    validModes = set(["alignment", "deepTools_qc"])
    for mode in modeTemp:
        if mode not in validModes:
            sys.exit("{} is not a valid mode!\n".format(mode))
    if "alignment" not in modeTemp and args.UMIDedup:
        sys.exit("UMIDedup is only valid for \"alignment\" mode!\n")
    if args.fromBAM:
        args.aligner = "EXTERNAL_BAM"
    ## End workflow-specific checks

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.plot_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
