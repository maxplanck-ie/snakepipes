import argparse
import os.path
import glob


def ListGenomes():
    """
    Return a list of all genome yaml files (sans the .yaml suffix)
    """
    dName = os.path.dirname(__file__)
    genomes = [os.path.basename(f)[:-5] for f in glob.glob(os.path.join(dName, "shared/organisms/*.yaml"))]
    return genomes
    
    
def mainArguments(defaults, workingDir=False):
    """
    Return a parser with the general and required args. This will include EITHER
    a -d option OR -i and -o, depending on the workingDir setting

    defaults is a dictionary of default values
    """

    parser = argparse.ArgumentParser(add_help=False)

    genomes = ListGenomes()
    parser.add_argument("genome", metavar="GENOME", help="Genome acronym of the target organism. Either a yaml file or one of: {}".format(", ".join(genomes)))

    required = parser.add_argument_group('required arguments')
    if workingDir:
        required.add_argument("-d", "--working-dir",
                              dest="workingdir",
                              help="working directory is output directory and must contain DNA-mapping pipeline output files",
                              required=True)
    else:
        required.add_argument("-i", "--input-dir",
                              dest="indir",
                              required=True,
                              help="input directory containing the FASTQ files, either paired-end OR single-end data")
        required.add_argument("-o", "--output-dir",
                              dest="outdir",
                              required=True,
                              help="output directory")
    

    general = parser.add_argument_group('general arguments')
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

    general.add_argument("--snakemake_options",
                         dest="snakemake_options",
                         metavar="STR",
                         type=str,
                         help="Snakemake options to be passed directly to snakemake, e.g. use --snakemake_options='--dryrun --rerun-incomplete --unlock --forceall'. (default: '%(default)s')",
                         default=defaults["snakemake_options"])

    general.add_argument("--tempdir",
                         dest="tempdir",
                         type=str,
                         help="used prefix path for temporary directory created via mktemp. Created temp dir gets exported as $TMPDIR and is removed at the end of this wrapper! (default: '%(default)s')",
                         default=defaults["tempdir"])

    general.add_argument("--DAG",
                         dest="createDAG",
                         action="store_true",
                         help="If specified, a file called pipeline.pdf is produced in the output directory that shows the rules used and their relationship to each other.")

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
