__description__ = """
MPI-IE workflow for ChIPseq analysis

Usage example:
    ChIPseq -d working-dir mm10 samples.yaml
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
                         "pairedEnd": True, "fragmentLength":200, "bwBinSize": 25,
                         "windowSize": 150, "predictChIPDict": None, "fromBAM": False,
                         "bigWigType": "both", "peakCaller": "MACS2", "sampleSheet": "",
                         "externalBed": "",
                         "plotFormat": "png", "bamExt": ".filtered.bam", "fdr": 0.05,
                         "absBestLFC": 1, "useSpikeinForNorm": False, "spikeinExt": "_spikein",
                         "peakCallerOptions": "--qvalue 0.001","cutntag": False,
                         "getSizeFactorsFrom": "genome"}):

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

    parser.add_argument("samples_config", nargs='?', metavar="SAMPLESCONFIG", help="configuration file (eg. 'example.chip_samples.yaml') with sample annotation")


    # Workflow options
    optional = parser.add_argument_group('Options')
    optional.add_argument("--peakCaller",
                          help="The peak caller to use. The default is %(default)s and this is only applicable for sharper peaks (broad peaks will always use histoneHMM).",
                          choices=["MACS2", "Genrich", "SEACR"],
                          default=defaults['peakCaller'])

    optional.add_argument("--peakCallerOptions",
                          help="Custom options to set for the peak caller. Default is '%(default)s'.",
                          default=defaults['peakCallerOptions'])

    optional.add_argument("--cutntag",
                          help="if set, MACS2 peakCaller is used with the parameters have been used "
                          "in the method section of Meers et al. 2019, Kaya-Okur et al. 2019 and 2020. "
                          "Setting this flag overwrites the '--peakCallerOptions'."
                          " Default is '%(default)s'.",
                          action="store_true")

    optional.add_argument("--singleEnd",
                          dest="pairedEnd",
                          action="store_false",
                          help="Input data is single-end, not paired-end")

    optional.add_argument("--useSpikeInForNorm",
                          dest="useSpikeInForNorm",
                          action="store_true",
                          help="Use the spikeIn chromosomes of the hybrid genome for normalization.")

    optional.add_argument("--getSizeFactorsFrom",
                          dest="getSizeFactorsFrom",
                          action="store",
                          choices=["genome", "TSS", "input"],
                          help="Which part of the spikein genome to use to calculate sizeFactors from.",
                          default=defaults['getSizeFactorsFrom'])

    optional.add_argument("--spikeinExt",
                          dest="spikeinExt",
                          help="Extention of spikein chromosome names in the hybrid genome. Ignored if useSpikeInForNorm is False (default: '%(default)s') .",
                          default=defaults["spikeinExt"])

    optional.add_argument("--bigWigType",
                          help="Type of bigWig file to create. Options are: 'subtract' (control-subtracted ChIP coverage), 'log2ratio' (for log2 ratio of ChIP over control) or 'both' (create both set of bed files). Note that the allele-specific mode currently only produces 'log2ratio' bigwigs. (default: '%(default)s')",
                          default=defaults["bigWigType"])

    optional.add_argument("--fragmentLength",
                          metavar="INT",
                          type=int,
                          help="Fragment length in sequencing. Used only if --singleEnd."
                               "(default: '%(default)s')",
                          default=defaults["fragmentLength"])

    optional.add_argument("--bwBinSize",
                          metavar="INT",
                          help="bin size of output files in bigWig format (default: '%(default)s')",
                          type=int,
                          default=defaults["bwBinSize"])

    optional.add_argument("--sampleSheet",
                          help="Information on samples (If differential binding analysis required); see "
                               "'https://github.com/maxplanck-ie/snakepipes/tree/master/docs/content/sampleSheet.example.tsv' for example. "
                               "IMPORTANT: The first entry defines which group of samples are control. "
                               "By this, the order of comparison and likewise the sign of values can be changed! "
                               "Also, the condition `control` should only be used for input samples (control peaks "
                               "are not evaluated for differential binding) (default: '%(default)s')",
                          default=defaults["sampleSheet"])

    optional.add_argument("--externalBed",
                          help="A bed file with intervals to be tested for differential binding. (default: '%(default)s')",
                          default=defaults["externalBed"])

    optional.add_argument("--windowSize",
                          help="Window size to counts reads in (If differential binding analysis required); "
                               "Default size is suitable for most transcription factors and sharp histone marks. "
                               "Small window sizes (~20bp) should be used for very narrow transcription factor peaks, "
                               "while large window sizes (~500 bp) should be used for broad marks (eg. H3K27me3) "
                               "(default: '%(default)s')",
                          default=defaults["windowSize"])

    optional.add_argument("--predictChIPDict",
                          nargs='?',
                          action='store',
                          help="Use existing bam files to predict a ChIPseq sample configuration file. Write it to the workingdir. "
                               "If no value is given, samples that contain 'input' are used as ChIP input/ctrl. Provide a custom pattern like 'input,H3$,H4$' to change that!",
                          default= None,
                          const="input")

    optional.add_argument("--fromBAM",
                          dest="fromBAM",
                          help="Input folder with bam files. If provided, the analysis will start from this point. If bam files contain single ends, please specify --singleEnd additionally. (default: '%(default)s')",
                          default=defaults["fromBAM"])

    optional.add_argument("--bamExt",
                          help="Extention of provided bam files, will be substracted from basenames to obtain sample names. (default: '%(default)s')",
                          default=defaults["bamExt"])


    optional.add_argument("--plotFormat",
                         choices=['png', 'pdf', 'None'],
                         help="Format of the output plots from deepTools. Select 'none' for no plots (default: '%(default)s')",
                         default=defaults["plotFormat"])

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

    # Local argument checks
    if args.predictChIPDict is not None:
        cf.predict_chip_dict(args.workingdir, args.predictChIPDict, args.bamExt, args.fromBAM)
        sys.exit(0)

    if args.samples_config is not None and os.path.exists(os.path.abspath(args.samples_config)):
        args.samples_config = os.path.abspath(args.samples_config)
    else:
        sys.exit("\nError! Sample config file not found! ({})\n Either create one or call the workflow with --predictChIPDict!\n".format(args.samples_config))

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.print_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
