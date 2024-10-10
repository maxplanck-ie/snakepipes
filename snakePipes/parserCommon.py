import argparse
import os.path
import glob
from importlib.metadata import version


def ListGenomes():
    """
    Return a list of all genome yaml files (sans the .yaml suffix)
    """
    dName = os.path.dirname(__file__)
    genomes = [os.path.basename(f)[:-5] for f in glob.glob(os.path.join(dName, "shared/organisms/*.yaml"))]
    return genomes


def mainArguments(defaults, workingDir=False, createIndices=False, preprocessing=False):
    """
    Return a parser with the general and required args. This will include EITHER
    a -d option OR -i and -o, depending on the workingDir setting

    defaults is a dictionary of default values

    A number of standard arguments are eliminated in the createIndices workflow.
    """

    # Set up some defaults for the sake of readthedocs
    if 'smtpServer' not in defaults:
        defaults['smtpServer'] = None
    if 'smtpPort' not in defaults:
        defaults['smtpPort'] = 0
    if 'onlySSL' not in defaults:
        defaults['onlySSL'] = False
    if 'emailSender' not in defaults:
        defaults['emailSender'] = None

    parser = argparse.ArgumentParser(add_help=False)

    if not createIndices and not preprocessing:
        genomes = ListGenomes()
        parser.add_argument("genome", metavar="GENOME", help="Genome acronym of the target organism. Either a yaml file or one of: {}".format(", ".join(genomes)))

    required = parser.add_argument_group('Required Arguments')
    if workingDir:
        required.add_argument("-d", "--working-dir",
                              dest="workingdir",
                              help="working directory is output directory and must contain DNAmapping pipeline output files",
                              required=True)
    else:
        if not createIndices:
            required.add_argument("-i", "--input-dir",
                                  dest="indir",
                                  required=True,
                                  help="input directory containing the FASTQ files, either paired-end OR single-end data")
        required.add_argument("-o", "--output-dir",
                              dest="outdir",
                              required=True,
                              help="output directory")

    general = parser.add_argument_group('General Arguments')
    general.add_argument("-h", "--help",
                         action="help",
                         help="show this help message and exit")

    general.add_argument("-v", "--verbose",
                         dest="verbose",
                         action="store_true",
                         help="verbose output (default: '%(default)s')",
                         default=defaults["verbose"])

    if not workingDir and not createIndices:
        general.add_argument("--ext",
                             help="Suffix used by input fastq files (default: '%(default)s').",
                             default=defaults["ext"])

        general.add_argument("--reads",
                             nargs=2,
                             help="Suffix used to denote reads 1 and 2 for paired-end data. This should typically be either '_1' '_2' or '_R1' '_R2' (default: '%(default)s). "
                             "Note that you should NOT separate the values by a comma (use a space) or enclose them in brackets.",
                             default=defaults["reads"])

    general.add_argument("-c", "--configFile",
                         help="configuration file: config.yaml (default: '%(default)s')",
                         default=defaults["configFile"])

    general.add_argument("--keepTemp",
                         action="store_true",
                         help="Prevent snakemake from removing files marked as being temporary (typically intermediate files that are rarely needed by end users). This is mostly useful for debugging problems.")

    general.add_argument("--snakemakeOptions",
                         action="append",
                         help="Snakemake options to be passed directly to snakemake, e.g. use --snakemakeOptions='--dryrun --rerun-incomplete --unlock --forceall'. WARNING! ONLY EXPERT USERS SHOULD CHANGE THIS! THE DEFAULT VALUE WILL BE APPENDED RATHER THAN OVERWRITTEN! (default: '%(default)s')",
                         default=[defaults["snakemakeOptions"]])

    general.add_argument("--DAG",
                         dest="createDAG",
                         action="store_true",
                         help="If specified, a file ending in _pipeline.pdf is produced in the output directory that shows the rules used and their relationship to each other.")

    general.add_argument("--version",
                         action="version",
                         version="%(prog)s {}".format(version("snakePipes")))

    emailArgs = parser.add_argument_group('Email Arguments')
    emailArgs.add_argument("--emailAddress",
                           help="If specified, send an email upon completion to the given email address")

    emailArgs.add_argument("--smtpServer",
                           default=defaults["smtpServer"],
                           help="If specified, the email server to use.")

    emailArgs.add_argument("--smtpPort",
                           type=int,
                           default=defaults["smtpPort"],
                           help="The port on the SMTP server to connect to. A value of 0 specifies the default port.")

    emailArgs.add_argument("--onlySSL",
                           action="store_true",
                           default=defaults["onlySSL"],
                           help="The SMTP server requires an SSL connection from the beginning.")

    emailArgs.add_argument("--emailSender",
                           default=defaults["emailSender"],
                           help="The address of the email sender. If not specified, it will be the address indicated by `--emailAddress`")

    emailArgs.add_argument("--smtpUsername",
                           help="If your SMTP server requires authentication, this is the username to use.")

    emailArgs.add_argument("--smtpPassword",
                           help="If your SMTP server requires authentication, this is the password to use.")

    return parser


def snpArguments(defaults):
    """
    Arguments related to allele-specific pipelines
    """
    parser = argparse.ArgumentParser(add_help=False)
    snpargs = parser.add_argument_group('Allele-specific mapping arguments')
    snpargs.add_argument("--VCFfile",
                         default='',
                         help="VCF file to create N-masked genomes (default: 'None'). Note that for the makePairs workflow this file is assumed to be gzipped and indexed (with tabix).")

    snpargs.add_argument("--strains",
                         default='',
                         help="Name or ID of SNP strains separated by comma (default: 'None')")

    snpargs.add_argument("--SNPfile",
                         default='',
                         help="File containing SNP locations (default: 'None')")

    snpargs.add_argument("--NMaskedIndex",
                         default='',
                         help="N-masked index of the reference genome (default: 'None'). "
                         "Note that this should point to a file (i.e. 'Genome' for STAR indices, genome.1.bt2 for bowtie2 indices).")

    return parser


# DNAmapping options added
def commonOptions(grp, defaults, bw=True, plots=True, preprocessing=False):
    """
    Common options found in many workflows
    grp is an argument group that's simply appended to
    """

    if not preprocessing:
        grp.add_argument("--downsample",
                         dest="downsample",
                         metavar="INT",
                         help="Downsample the given number of reads randomly from of each FASTQ file (default: '%(default)s')",
                         type=int,
                         default=defaults["downsample"])

        grp.add_argument("--trim",
                         dest="trim",
                         action="store_true",
                         help="Activate fastq read trimming. If activated, Illumina adaptors are trimmed by default. "
                         "Additional parameters can be specified under --trimmerOptions. (default: '%(default)s')",
                         default=defaults["trim"])

        grp.add_argument("--trimmer",
                         dest="trimmer",
                         choices=['cutadapt', 'trimgalore', 'fastp'],
                         help="Trimming program to use: Cutadapt, TrimGalore, or fastp. Note that if you change this you may "
                         "need to change --trimmerOptions to match! (default: '%(default)s')",
                         default=defaults["trimmer"])

        grp.add_argument("--trimmerOptions",
                         dest="trimmerOptions",
                         help="Additional option string for trimming program of choice. (default: '%(default)s')",
                         default=defaults["trimmerOptions"])

    grp.add_argument("--fastqc",
                     dest="fastqc",
                     action="store_true",
                     help="Run FastQC read quality control (default: '%(default)s')",
                     default=defaults["fastqc"])

    grp.add_argument("--bcExtract",
                     dest="UMIBarcode",
                     action="store_true",
                     help="To extract umi barcode from fastq file via UMI-tools and add it to the read name "
                     "(default: '%(default)s')",
                     default=defaults["UMIBarcode"])

    grp.add_argument("--bcPattern",
                     help="The pattern to be considered for the barcode. 'N' = UMI position (required) 'C' = barcode position (optional) "
                     "(default: '%(default)s')",
                     default=defaults["bcPattern"])

    if not preprocessing:
        grp.add_argument("--UMIDedup",
                         action="store_true",
                         help="Deduplicate bam file based on UMIs via `umi_tools dedup` that are present in the read name. "
                         "(default: '%(default)s')",
                         default=defaults["UMIDedup"])

        grp.add_argument("--UMIDedupSep",
                         help="umi separation character "
                         "that will be passed to umi_tools."
                         "(default: '%(default)s')",
                         default=defaults["UMIDedupSep"])

        grp.add_argument("--UMIDedupOpts",
                         help="Additional options that will be passed to umi_tools."
                         "(default: '%(default)s')",
                         default=defaults["UMIDedupOpts"])

    if bw and not preprocessing:
        grp.add_argument("--bwBinSize",
                         dest="bwBinSize",
                         help="Bin size of output files in bigWig format (default: '%(default)s')",
                         type=int,
                         default=defaults["bwBinSize"])

    if plots and not preprocessing:
        grp.add_argument("--plotFormat",
                         choices=['png', 'pdf', 'None'],
                         metavar="STR",
                         type=str,
                         help="Format of the output plots from deepTools. Select 'none' for no plots (default: '%(default)s')",
                         default=defaults["plotFormat"])
