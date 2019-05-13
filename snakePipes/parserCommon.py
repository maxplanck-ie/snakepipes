import argparse
import os.path
import glob
from snakePipes import __version__


def ListGenomes():
    """
    Return a list of all genome yaml files (sans the .yaml suffix)
    """
    dName = os.path.dirname(__file__)
    genomes = [os.path.basename(f)[:-5] for f in glob.glob(os.path.join(dName, "shared/organisms/*.yaml"))]
    return genomes


def mainArguments(defaults, workingDir=False, createIndices=False):
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

    if not createIndices:
        genomes = ListGenomes()
        parser.add_argument("genome", metavar="GENOME", help="Genome acronym of the target organism. Either a yaml file or one of: {}".format(", ".join(genomes)))

    required = parser.add_argument_group('Required Arguments')
    if workingDir:
        required.add_argument("-d", "--working-dir",
                              dest="workingdir",
                              help="working directory is output directory and must contain DNA-mapping pipeline output files",
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

    general.add_argument("-c", "--configfile",
                         dest="configfile",
                         help="configuration file: config.yaml (default: '%(default)s')",
                         default=defaults["configfile"])

    general.add_argument("--cluster_configfile",
                         dest="cluster_configfile",
                         help="configuration file for cluster usage. In absence, the default options "
                         "from shared/cluster.yaml and workflows/[workflow]/cluster.yaml would be selected (default: '%(default)s')",
                         default=defaults["cluster_configfile"])

    general.add_argument("-j", "--jobs",
                         dest="max_jobs",
                         metavar="INT",
                         help="maximum number of concurrently submitted Slurm jobs / cores if workflow is run locally (default: '%(default)s')",
                         type=int, default=defaults["max_jobs"])

    general.add_argument("--local",
                         dest="local",
                         action="store_true",
                         default=False,
                         help="run workflow locally; default: jobs are submitted to Slurm queue (default: '%(default)s')")

    general.add_argument("--keepTemp",
                         action="store_true",
                         help="Prevent snakemake from removing files marked as being temporary (typically intermediate files that are rarely needed by end users). This is mostly useful for debugging problems.")

    general.add_argument("--snakemake_options",
                         dest="snakemake_options",
                         metavar="STR",
                         type=str,
                         action="append",
                         help="Snakemake options to be passed directly to snakemake, e.g. use --snakemake_options='--dryrun --rerun-incomplete --unlock --forceall'. WARNING! ONLY EXPERT USERS SHOULD CHANGE THIS! THE DEFAULT VALUE WILL BE APPENDED RATHER THAN OVERWRITTEN! (default: '%(default)s')",
                         default=[defaults["snakemake_options"]])

    general.add_argument("--tempdir",
                         dest="tempdir",
                         type=str,
                         help="used prefix path for temporary directory created via mktemp. Created temp dir gets exported as $TMPDIR and is removed at the end of this wrapper! (default: '%(default)s')",
                         default=defaults["tempdir"])

    general.add_argument("--DAG",
                         dest="createDAG",
                         action="store_true",
                         help="If specified, a file ending in _pipeline.pdf is produced in the output directory that shows the rules used and their relationship to each other.")

    general.add_argument("--version",
                         action="version",
                         version="%(prog)s {}".format(__version__))

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
                         metavar="STR",
                         help="VCF file to create N-masked genomes (default: 'None')",
                         default='',
                         type=str)

    snpargs.add_argument("--strains",
                         dest="strains",
                         metavar="STR",
                         help="Name or ID of SNP strains separated by comma (default: 'None')",
                         default='',
                         type=str)

    snpargs.add_argument("--SNPfile",
                         metavar="STR",
                         help="File containing SNP locations (default: 'None')",
                         default='',
                         type=str)

    snpargs.add_argument("--Nmasked_index",
                         metavar="STR",
                         help="N-masked index of the reference genome (default: 'None')",
                         default='',
                         type=str)

    return parser


# DNA-mapping options added
def commonOptions(grp, defaults, bw=True, plots=True):
    """
    Common options found in many workflows
    grp is an argument group that's simply appended to
    """
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
                     "Additional parameters can be specified under --trim_options. (default: '%(default)s')",
                     default=defaults["trim"])

    grp.add_argument("--trim_prg",
                     dest="trim_prg",
                     choices=['cutadapt', 'trimgalore'],
                     help="Trimming program to use: Cutadapt or TrimGalore. Note that if you change this you may "
                     "need to change --trim_options to match! (default: '%(default)s')",
                     default=defaults["trim_prg"])

    grp.add_argument("--trim_options",
                     dest="trim_options",
                     metavar="STR",
                     type=str,
                     help="Additional option string for trimming program of choice. (default: '%(default)s')",
                     default=defaults["trim_options"])

    grp.add_argument("--fastqc",
                     dest="fastqc",
                     action="store_true",
                     help="Run FastQC read quality control (default: '%(default)s')",
                     default=defaults["fastqc"])

    if bw:
        grp.add_argument("--bw-binsize",
                         dest="bw_binsize",
                         metavar="INT",
                         help="Bin size of output files in bigWig format (default: '%(default)s')",
                         type=int,
                         default=defaults["bw_binsize"])

    if plots:
        grp.add_argument("--plotFormat",
                         dest="plot_format",
                         choices=['png', 'pdf', 'None'],
                         metavar="STR",
                         type=str,
                         help="Format of the output plots from deepTools. Select 'none' for no plots (default: '%(default)s')",
                         default=defaults["plot_format"])
