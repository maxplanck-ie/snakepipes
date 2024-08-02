#!/usr/bin/env python3

__description__ = """
MPI-IE workflow for RNA mapping and analysis

usage example:
    RNAseq -i input-dir -o output-dir mm10
"""


import argparse
import os
import sys
import textwrap
import snakePipes.common_functions as cf
import snakePipes.parserCommon as parserCommon
import warnings


def parse_args(defaults={"verbose": False, "configFile": None,
                         "clusterConfigFile": None, "maxJobs": 5,
                         "snakemakeOptions": "--use-conda", "tempDir": None,
                         "mode": "alignment,deepTools_qc", "downsample": False,
                         "trim": False, "trimmer": "cutadapt",
                         "trimmerOptions": None, "fastqc": False,
                         "libraryType": 2, "aligner": "STAR",
                         "alignerOptions": None,
                         "featureCountsOptions": "-C -Q 10 --primary",
                         "filterGTF": None, "sampleSheet": None,
                         "formula": "",
                         "reads": ["_R1", "_R2"], "ext": ".fastq.gz",
                         "bwBinSize": 25, "dnaContam": False, "plotFormat": "png",
                         "fromBAM": False, "bamExt": ".bam", "pairedEnd": True,
                         "UMIDedup": False,
                         "UMIDedupOpts": "", "bcPattern": "NNNNCCCCCCCCC",
                         "UMIDedupSep": "_", "UMIBarcode": False, "rMats": False,
                         "fdr": 0.05}):
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
                          help="workflow running modes (available: 'alignment-free, alignment, allelic-mapping, allelic-counting, deepTools_qc, three-prime-seq')"
                          " (default: '%(default)s')",
                          default=defaults["mode"])

    parserCommon.commonOptions(optional, defaults, bw=True)

    optional.add_argument("--libraryType",
                          help="user provided library type strand specificity. featureCounts style: 0, 1, 2 (Illumina TruSeq); default: '%(default)s')",
                          type=int,
                          default=defaults["libraryType"])

    optional.add_argument("--aligner",
                          help="Program used for mapping: STAR or HISAT2 (default: '%(default)s'). If you change this, please change --alignerOptions to match.",
                          choices=["STAR","HISAT2"],
                          default=defaults["aligner"])

    optional.add_argument("--alignerOptions",
                          help="STAR or hisat2 option string, e.g.: '--twopassMode Basic' (default: '%(default)s')",
                          default=defaults["alignerOptions"])


    optional.add_argument("--featureCountsOptions",
                          help="featureCounts option string. The options '-p -B'"
                               " are always used for paired-end data (default: '%(default)s')",
                          default=defaults["featureCountsOptions"])

    optional.add_argument("--filterGTF",
                          help="filter annotation GTF by grep for use with Salmon, e.g."
                               " use --filterGTF='-v pseudogene'; default: '%(default)s')",
                          default=defaults["filterGTF"])

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


    optional.add_argument("--dnaContam",
                          action="store_true",
                          help="Returns a plot which presents the proportion of the intergenic reads (default: '%(default)s')",
                          default=defaults["dnaContam"])

    optional.add_argument("--fromBAM",
                         action="store_true",
                         help="Input folder with bam files. If provided, the analysis will start from this point. If bam files contain single ends, please specify --singleEnd additionally.",
                         default=defaults["fromBAM"])

    optional.add_argument("--bamExt",
                          help=("Extention of provided bam files, will be "
                          "subtracted from basenames to obtain sample names. "
                          "(default: '%(default)s')"),
                          default=defaults["bamExt"])


    optional.add_argument("--singleEnd",
                          dest="pairedEnd",
                          action="store_false",
                          help="input data is single-end, not paired-end. This is only used if --fromBAM is specified.")

    optional.add_argument("--rMats",
                          dest="rMats",
                          action="store_true",
                          help="Run differential splicing analysis using rMats-turbo. Note that this flag requires --sampleSheet to be specified.",
                          default=defaults["rMats"])

    optional.add_argument("--FDR",
                          dest="fdr",
                          help="FDR threshold to apply for filtering DE genes"
                               "(default: '%(default)s')",
                          default=defaults["fdr"])

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
    # check for Allele-specific mapping mode
    args.allele_mode = cf.checkAlleleParams(args)
    # convert file path to abspath
    if args.allele_mode == "create_and_map":
        args.VCFfile = os.path.abspath(args.VCFfile)
    elif args.allele_mode == "map_only":
        args.SNPfile = os.path.abspath(args.SNPfile)
        args.NMaskedIndex = os.path.abspath(args.NMaskedIndex)
    modeTemp = args.mode.split(",")
    validModes = set(["alignment", "alignment-free", "deepTools_qc", "allelic-mapping", "allelic-counting", "three-prime-seq"])
    for mode in modeTemp:
        if mode not in validModes:
            sys.exit("{} is not a valid mode!\n".format(mode))
    if "alignment" not in modeTemp and args.UMIDedup:
        sys.exit("UMIDedup is only valid for \"alignment\" mode!\n")
    if "allelic-counting" in modeTemp and "deepTools_qc" in modeTemp:
        sys.exit("Mode deepTools_qc is not compatible with mode allelic-counting.")
    if args.fromBAM and ("alignment-free" in modeTemp ):
        sys.exit("\n--fromBAM can only be used with modes \'alignment\' , \'allelic-mapping\' , \'allelic-counting\'  or \'deepTools_qc\' - use one of these modes or provide fastq files!\n")
    if args.fromBAM:
        args.aligner = "EXTERNAL_BAM"
    if "allelic-counting" in modeTemp and not args.fromBAM:
        warnings.warn("--fromBAM is required with allelic-counting mode. Setting to True.")
        args.fromBAM = True
    if "allelic-counting" in modeTemp:
        args.bamExt = ".sorted.bam"
        args.aligner = "allelic_bams"
    if args.rMats and not args.sampleSheet:
        sys.exit("--rMats flag requires a sampleSheet (specified with --sampleSheet).\n")
    if "three_prime_seq" in mode:
        if not args.sampleSheet:
            sys.exit("mode three-prime-seq requires a sampleSheet "
                     "(specified with --sampleSheet).\n")
        aligner = "STAR"
        alignerOptions = defaults['threePrimeAlignerOptions']
        trimmerOptions = defaults['threePrimeTrimmerOptions']
        trimmer = "fastp"
        trim = True
    
    ## End workflow-specific checks

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.print_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
