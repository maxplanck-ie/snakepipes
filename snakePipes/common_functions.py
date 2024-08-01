#!/usr/bin/env python

# functions shared across workflows ##########################################
##############################################################################
import subprocess
import os
import re
import yaml
import glob
import sys
import shutil
from pathlib import Path
from thefuzz import fuzz
import smtplib
from email.message import EmailMessage
#from snakePipes import __version__
from importlib.metadata import version


def set_env_yamls():
    """
    This defines the global variables describing where the conda env yaml files are
    """
    return {'CONDA_SHARED_ENV': 'envs/shared.yaml',
            'CONDA_CREATE_INDEX_ENV': 'envs/createIndices.yaml',
            'CONDA_RNASEQ_ENV': 'envs/rna_seq.yaml',
            'CONDA_SALMON_ENV': 'envs/salmon.yaml',
            'CONDA_SLEUTH_ENV': 'envs/sleuth.yaml',
            'CONDA_RMATS_ENV': 'envs/rMats.yaml',
            'CONDA_scRNASEQ_ENV': 'envs/sc_rna_seq.yaml',
            'CONDA_seurat_ENV': 'envs/sc_rna_seq_seurat.yaml',
            'CONDA_loompy_ENV': 'envs/sc_rna_seq_loompy.yaml',
            'CONDA_alevinqc_ENV': 'envs/sc_rna_seq_alevinqc.yaml',
            'CONDA_eisaR_ENV': 'envs/sc_rna_seq_eisaR.yaml',
            'CONDA_DNA_MAPPING_ENV': 'envs/dna_mapping.yaml',
            'CONDA_CHIPSEQ_ENV': 'envs/chip_seq.yaml',
            'CONDA_ATAC_ENV': 'envs/atac_seq.yaml',
            'CONDA_HIC_ENV': 'envs/hic.yaml',
            'CONDA_MAKEPAIRS_ENV': 'envs/makePairs.yaml',
            'CONDA_WGBS_ENV': 'envs/wgbs.yaml',
            'CONDA_DSS_ENV': 'envs/wgbs_dss.yaml',
            'CONDA_RMD_ENV': 'envs/rmarkdown.yaml',
            'CONDA_PREPROCESSING_ENV': 'envs/preprocessing.yaml',
            'CONDA_NONCODING_RNASEQ_ENV': 'envs/noncoding.yaml',
            'CONDA_SAMBAMBA_ENV': 'envs/sambamba.yaml',
            'CONDA_pysam_ENV': 'envs/pysam.yaml',
            'CONDA_SEACR_ENV': 'envs/chip_seacr.yaml'}


def merge_dicts(x, y):
    z = {}
    z = x.copy()
    if y:
        z.update(y)
    return z


# this is a pure sanity function to avoid obvious mailfunction during snakefile execution
# because we load yaml/path/genome configs directly into global namespace!
def sanity_dict_clean(myDict):
    unwanted_keys = ['maindir', 'workflow']
    for k in unwanted_keys:
        if myDict and k in myDict:
            del myDict[k]
    return myDict


def namesOKinR(sampleNames):
    """
    Return nothing, but print warning to the screen
    if any of the sample names will get munged by R.
    """
    reservedWords = set(["NULL", "NA", "TRUE", "FALSE", "Inf", "NaN", "NA_integer_", "NA_real_",
                         "NA_character_", "NA_complex_", "function", "while", "repeat", "for",
                         "if", "in", "else", "next", "break", "..."])
    for sampleName in sampleNames:
        # Starts with A-Za-z or .
        if (not sampleName[0].isalpha()) and (not sampleName[0] == "."):
            sys.stderr.write("Any steps involving R packages will fail if sample names do not start with a letter or '.'. {} is not compatible and will fail these!\n".format(sampleName))
        # reserved word
        if sampleName in reservedWords:
            sys.stderr.write("{} is a reserved keyword in R, so if there are steps using R they will fail!\n".format(sampleName))
        # invalid characters, which is everything except alpha numeric, . and _
        if not all([(x.isalnum() or x in ["_", "."]) for x in sampleName]):
            sys.stderr.write("R requires that all samples names contain ONLY letters, number, '_' or '.', so {} is invalid and may cause failure in steps using R!\n".format(sampleName))


def load_configfile(configFiles, verbose, info='Config'):
    with open(configFiles, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = sanity_dict_clean(config)

    if verbose:
        print("\n--- " + info + " ---------------------------------------------------------------------")
        print("config file: {}".format(configFiles))
        for k, v in sorted(config.items()):
            print("{}: {}".format(k, v))
        print("-" * 80, "\n")
    return config


def write_configfile(configFile, config):
    with open(configFile, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)


# returns all key-value pairs that are different from dict1 to dict2
def config_diff(dict1, dict2):
    diff = {}
    for k in dict1:
        if k in dict2:
            if dict1[k] != dict2[k]:
                diff[k] = dict1[k]
        else:
            diff[k] = dict1[k]
    return diff


def get_version():
    # If this is sent to stdout it breaks making a DAG pdf
    sys.stderr.write("\n---- This analysis has been done using snakePipes version {} ----\n".format(version("snakePipes")))


def load_organism_data(genome, maindir, verbose):
    # Load the global config file, which dictates where the organisms should be found
    cfg = load_configfile(os.path.join(maindir, "shared", "defaults.yaml"), False, "defaults")

    if os.path.isfile(os.path.join(maindir, cfg['organismsDir'], genome + ".yaml")):
        organism = load_configfile(os.path.join(maindir, cfg['organismsDir'], genome + ".yaml"), verbose, "Genome")
    elif os.path.isfile(os.path.join(cfg['organismsDir'], genome + ".yaml")):
        organism = load_configfile(os.path.join(cfg['organismsDir'], genome + ".yaml"), verbose, "Genome")
    elif os.path.isfile(genome):
        organism = load_configfile(genome, verbose, "Genome (user)")
    else:
        exit("ERROR: Genome configuration file NOT found for: {}\n".format(genome))
    return organism


def get_sample_names(infiles, ext, reads):
    """
    Get sample names without file extensions
    """
    s = set()
    lext = len(ext)
    l0 = len(reads[0])
    l1 = len(reads[1])
    for x in infiles:
        x = os.path.basename(x)[:-lext]
        if x.endswith(reads[0]):
            x = x[:-l0]
            s.add(x)
        elif x.endswith(reads[1]):
            x = x[:-l1]
            s.add(x)
        else:
            sys.stderr.write("Warning! {} does not have {} as its name suffix. "
                             "Either change it or modify the 'reads' in the "
                             "config.yaml to your deired ones.\n".format(x, reads))

    if sorted(list(s)) == []:
        sys.exit("Error! No sample has the right read suffix ({}). "
                 "Please modify them or update the config.yaml with "
                 "your desired suffix.".format(reads))
    return sorted(list(s))


def get_sample_names_bam(infiles, bamExt):
    """
    Get sample names without file extensions
    """
    s = []
    for x in infiles:
        x = os.path.basename(x).replace(bamExt, "")
        s.append(x)
    return sorted(list(set(s)))


def get_sample_names_suffix_bam(infiles, bamExt):
    """
    Get sample names without file extensions
    """
    bamSuff = [x + bamExt for x in [".genome1", ".genome2", ".unassigned", ".allele_flagged"]]
    s = []
    for x in infiles:
        for y in bamSuff:
            if y in os.path.basename(x):
                x = os.path.basename(x).replace(y, "")
                s.append(x)
    return sorted(list(set(s)))


def is_paired(infiles, ext, reads):
    """
    Check for paired-end input files
    """
    pairedEnd = False
    infiles_dic = {}
    for infile in infiles:
        fname = os.path.basename(infile).replace(ext, "")
        m = re.match("^(.+)(" + reads[0] + "|" + reads[1] + ")$", fname)
        if m:
            bname = m.group(1)
            if bname not in infiles_dic:
                infiles_dic[bname] = [infile]
            else:
                infiles_dic[bname].append(infile)
    if not infiles_dic:
        sys.exit("Error: No fastq file has been found to be checked.")
    values_length = [len(x) for x in infiles_dic.values()]
    if min(values_length) == 2:
        pairedEnd = True
    elif min(values_length) == 1 and max(values_length) == 2:
        sys.exit("Error: The directory contains a mixture of paired-end and single-end data!")
    return pairedEnd


def check_replicates(sample_info_file):
    """
    return True if each condition has at least 2 replicates
    this check is eg. necessary for sleuth
    """
    f = open(sample_info_file)
    conditionCol = None
    nCols = None
    d = dict()
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            if "condition" not in cols or "name" not in cols:
                sys.exit("ERROR: Please use 'name' and 'condition' as column headers in the sample info file ({})!\n".format(sample_info_file))
            conditionCol = cols.index("condition")
            nCols = len(cols)
            continue
        elif idx == 1:
            # Sometimes there's a column of row names, which lack a header
            if len(cols) != nCols and len(cols) - 1 != nCols:
                sys.exit("ERROR: there's a mismatch between the number of columns in the header and body of {}!\n".format(sample_info_file))
            if len(cols) - 1 == nCols:
                conditionCol += 1
        if not len(line.strip()) == 0:
            if cols[conditionCol] not in d:
                d[cols[conditionCol]] = 0
            d[cols[conditionCol]] += 1
    f.close()

    for k, v in d.items():
        if v < 2:
            sys.stderr.write("WARNING: The {} group has no replicates!\n".format(k))
            return False

    return True


def isMultipleComparison(sampleSheet):
    f = open(sampleSheet)
    nCols = None
    d = dict()
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            if "group" not in cols:
                return False
            comparisonGroupCol = cols.index("group")
            nCols = len(cols)
            continue
        elif idx == 1:
            # Sometimes there's a column of row names, which lack a header
            if len(cols) - 1 == nCols:
                comparisonGroupCol += 1
        if not len(line.strip()) == 0:
            if cols[comparisonGroupCol] not in d:
                d[cols[comparisonGroupCol]] = 0
            d[cols[comparisonGroupCol]] += 1
    f.close()

    if len(d) > 1:
        return True


def splitSampleSheet(sampleSheet, destination_pfx):
    f = open(sampleSheet)
    conditionCol = None
    nameCol = None
    comparisonGroupCol = None
    batchCol = None
    nCols = None
    d = dict()
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            conditionCol = cols.index("condition")
            nameCol = cols.index("name")
            comparisonGroupCol = cols.index("group")
            if "batch" in cols:
                batchCol = cols.index("batch")
            nCols = len(cols)
            continue
        elif idx == 1:
            # Sometimes there's a column of row names, which lack a header
            if len(cols) != nCols and len(cols) - 1 != nCols:
                sys.exit("ERROR: there's a mismatch between the number of columns in the header and body of {}!\n".format(sampleSheet))
            if len(cols) - 1 == nCols:
                conditionCol += 1
                nameCol += 1
                comparisonGroupCol += 1
                if batchCol:
                    batchCol += 1
            firstCondition = cols[conditionCol]
        if not len(line.strip()) == 0:
            if cols[comparisonGroupCol] not in d:
                d[cols[comparisonGroupCol]] = []
            if batchCol:
                d[cols[comparisonGroupCol]].append([cols[nameCol], cols[batchCol], cols[conditionCol]])
            else:
                d[cols[comparisonGroupCol]].append([cols[nameCol], cols[conditionCol]])

    f.close()

    if "All" in d.keys():
        if batchCol:
            allCondition = d["All"][0][2]
        else:
            allCondition = d["All"][0][1]
        if allCondition == firstCondition:
            d["All"].reverse()
    for k in d.keys():
        if k != "All" and "All" in d.keys():
            if allCondition == firstCondition:
                for x in d["All"]:
                    d[k].insert(0, x)
            else:
                d[k].extend(d['All'])

        outfile = os.path.join("splitSampleSheets", '.'.join([os.path.basename(destination_pfx), k, 'tsv']))
        with open(outfile, 'w') as of:
            if batchCol:
                of.write('name\tbatch\tcondition\n')
            else:
                of.write('name\tcondition\n')
            for item in d[k]:
                of.write('\t'.join(item) + '\n')

    return


def returnComparisonGroups(sampleSheet):
    f = open(sampleSheet)
    nCols = None
    d = dict()
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            if "group" not in cols:
                return False
            comparisonGroupCol = cols.index("group")
            nCols = len(cols)
            continue
        elif idx == 1:
            # Sometimes there's a column of row names, which lack a header
            if len(cols) - 1 == nCols:
                comparisonGroupCol += 1
        if not len(line.strip()) == 0:
            if cols[comparisonGroupCol] not in d:
                d[cols[comparisonGroupCol]] = 0
            d[cols[comparisonGroupCol]] += 1
    f.close()

    if "All" in d.keys():
        del d['All']

    return d.keys()


def sampleSheetGroups(sampleSheet, multipleComp):
    """
    Parse a sampleSheet and return a dictionary with keys the group and values the sample names
    """
    f = open(sampleSheet)
    conditionCol = None
    nameCol = None
    groupCol = None
    nCols = None
    d = dict()
    for idx, line in enumerate(f):
        cols = line.strip().split("\t")
        if idx == 0:
            if "condition" not in cols or "name" not in cols:
                sys.exit("ERROR: Please use 'name' and 'condition' as column headers in the sample info file ({})!\n".format(sampleSheet))
            conditionCol = cols.index("condition")
            nameCol = cols.index("name")
            if multipleComp:
                groupCol = cols.index("group")
            nCols = len(cols)
            continue
        elif idx == 1:
            # Sometimes there's a column of row names, which lack a header
            if len(cols) != nCols and len(cols) - 1 != nCols:
                sys.exit("ERROR: there's a mismatch between the number of columns in the header and body of {}!\n".format(sampleSheet))
            if len(cols) - 1 == nCols:
                conditionCol += 1
                nameCol += 1
                if multipleComp:
                    groupCol += 1
        if not len(line.strip()) == 0:
            if not multipleComp:
                if cols[conditionCol] not in d:
                    d[cols[conditionCol]] = []
                d[cols[conditionCol]].append(cols[nameCol])
            else:
                if cols[groupCol] not in d:
                    d[cols[groupCol]] = {}
                if cols[conditionCol] not in d[cols[groupCol]]:
                    d[cols[groupCol]][cols[conditionCol]] = []
                d[cols[groupCol]][cols[conditionCol]].append(cols[nameCol])
    if "All" in d.keys():
        for k in d.keys():
            if k not in "All":
                d[k][list(d["All"].keys())[0]] = []
                for x in d["All"].values():
                    # don't use append as this results in a list of lists and causes issues downstream
                    d[k][list(d["All"].keys())[0]] += x
        del d['All']
    f.close()
    return d


def make_temp_dir(tempDir, fallback_dir, verbose=False):
    try:
        output = subprocess.check_output("mktemp -d -p " + tempDir + "/ tmp.snakemake.XXXXXXXX", shell=True, stderr=subprocess.STDOUT)
        temp_path = output.decode().rstrip() + "/"
    except subprocess.CalledProcessError:
        try:
            print("\nFailed to create temp dir under temp path prefix (" + tempDir + ")! Try fallback: " + fallback_dir + "/ ...")
            output = subprocess.check_output("mktemp -d -p " + fallback_dir + "/ tmp.snakemake.XXXXXXXX", shell=True, stderr=subprocess.STDOUT)
            temp_path = output.decode().rstrip() + "/"
        except subprocess.CalledProcessError:
            print("\nAlso failed to create temp dir under fallback prefix (" + fallback_dir + "/)!")
            exit(1)
    if verbose:
        print("\ntemp dir created: " + temp_path)
    return temp_path


def checkAlleleParams(args):
    # first some sanity checks
    mode = list(map(str.strip, re.split(',|;', args.mode)))
    mode = [element.lower() for element in mode]
    if "allelic-mapping" in mode and "mapping" in mode:
        print("\nError! Please specify either allelic-mapping or mapping for option --mode! \n")
        exit(1)
    if "allelic-mapping" in args.mode:
        if not os.path.exists(args.SNPfile):
            # if no SNPfile, check for a VCF file
            if os.path.exists(args.VCFfile):
                # check for strain ID
                if args.strains == '':
                    print("\nError! Please specify strain ID to extract from given VCF file for Allele-specific mapping! ({})\n".format(args.VCFfile))
                    exit(1)
                else:
                    allele_mode = 'create_and_map'
            else:
                print("\nError! Please specify either VCF file or SNP file for Allele-specific mapping! \n")
                exit(1)
        # If SNP file is present, check whether genome index also exists
        elif not os.path.exists(os.path.dirname(args.NMaskedIndex)):
            print("\nError! Please specify an n-masked index file for Allele-specific mapping! \n")
            exit(1)
        else:
            allele_mode = 'map_only'
    else:
        allele_mode = None
    return allele_mode


def cleanLogs(d, cluster_config):
    """
    Remove all empty log files, both in cluster_logs/ and */logs/
    """
    if "snakePipes_cluster_logDir" in cluster_config:
        path = os.path.join(d, cluster_config["snakePipes_cluster_logDir"], "*")
        if re.search("^/", cluster_config["snakePipes_cluster_logDir"]):
            path = os.path.join(cluster_config["snakePipes_cluster_logDir"], "*")
        for f in glob.glob(path):
            s = os.stat(f)
            if s.st_size == 0:
                os.remove(f)
    for f in glob.glob(os.path.join(d, "*", "logs", "*")):
        s = os.stat(f)
        if s.st_size == 0:
            os.remove(f)


def check_sample_info_header(sampleSheet_file):
    """
    return True in case sample info file contains column names 'name' and 'condition'
    """
    if not os.path.isfile(sampleSheet_file):
        sys.exit("ERROR: Cannot find sample info file! (--sampleSheet {})\n".format(sampleSheet_file))
    sampleSheet_file = os.path.abspath(sampleSheet_file)
    ret = open(sampleSheet_file).read().split("\n")[0].split("\t")
    if "name" in ret and "condition" in ret:
        sys.stderr.write("Sample sheet found and header is ok!\n")
    else:
        sys.exit("ERROR: Please use 'name' and 'condition' as column headers in sample info file! Please use a tab-delimited file! ({})\n".format(sampleSheet_file))
    return sampleSheet_file


def setDefaults(fileName):
    """
    Set a number of variables used in the wrappers and the defaults
    """
    # Script-neutral paths
    baseDir = os.path.dirname(__file__)
    workflowDir = os.path.join(baseDir, "workflows", fileName)

    # defaults
    defaults = load_configfile(os.path.join(workflowDir, "defaults.yaml"), False)
    globalDefaults = load_configfile(os.path.join(baseDir, "shared/defaults.yaml"), False)
    defaults = merge_dicts(defaults, globalDefaults)
    return baseDir, workflowDir, defaults


def handleUserArgs(args, defaults, args_func):
    """
    If a user supplies a custom YAML file then that must then replace the defaults.
    However command line options need to take precedence, so update defaults and
    simply reparse the command line options (with args_func().parse_args())
    """
    if args.configFile:
        if not os.path.exists(args.configFile):
            sys.exit("\nError! Provided configFile (-c) not found! ({})\n".format(args.configFile))
        user_config = load_configfile(args.configFile, False)
        defaults = merge_dicts(defaults, user_config)
        parser = args_func(defaults)
        args = parser.parse_args()
    defaults.update(vars(args))
    return args, defaults


def sendEmail(args, returnCode):
    """
    Try to send an email to the user. Errors must be non-fatal.
    """
    try:
        msg = EmailMessage()
        msg['Subject'] = "Snakepipes completed"
        msg['From'] = args.emailSender
        msg['To'] = args.emailAddress
        if returnCode == 0:
            msg.set_content("The pipeline finished successfully\n")
        else:
            msg.set_content("The pipeline failed with exit code {}\n".format(returnCode))

        if args.onlySSL:
            s = smtplib.SMTP_SSL(args.smtpServer, port=args.smtpPort)
        else:
            s = smtplib.SMTP(args.smtpServer, port=args.smtpPort)
        if args.smtpUsername:
            s.login(args.smtpUsername, args.smtpPassword)
        s.send_message(msg)
        s.quit()
    except:
        sys.stderr.write("An error occured while sending the email.\n")
        pass


def checkCommonArguments(args, baseDir, outDir=False, createIndices=False, preprocessing=False):
    """
    Check the wrapper arguments

    The createIndices workflow disables some of this
    """
    # Some workflows use a working dir, others and outdir
    if outDir:
        args.outdir = os.path.abspath(args.outdir)
        args.workingdir = args.outdir
    else:
        args.workingdir = os.path.abspath(args.workingdir)

    # 1. Dir path
    if not createIndices:
        if outDir:
            if os.path.exists(args.indir):
                args.indir = os.path.abspath(args.indir)
            else:
                sys.exit("\nError! Input dir not found! ({})\n".format(args.indir))
        else:
            if "fromBAM" in args and args.fromBAM:
                if os.path.exists(args.fromBAM):
                    os.makedirs(args.workingdir, exist_ok=True)
                    args.workingdir = os.path.abspath(args.workingdir)
                    args.fromBAM = os.path.abspath(args.fromBAM)
                else:
                    sys.exit("\nError! Directory with bam files (--fromBAM) not found! ({})\n".format(args.fromBAM))
            else:
                if os.path.exists(args.workingdir):
                    args.workingdir = os.path.abspath(args.workingdir)
                else:
                    sys.exit("\nError! Working-dir (-d) dir not found! ({})\n".format(args.workingdir))
            args.outdir = args.workingdir
    # 2. Sample info file
    if 'sampleSheet' in args and args.sampleSheet and not preprocessing:
        args.sampleSheet = check_sample_info_header(args.sampleSheet)
    # 3. get abspath from user provided genome/organism file
    if not createIndices and not preprocessing:
        if not os.path.isfile(os.path.join(baseDir, "shared/organisms/{}.yaml".format(args.genome))) and os.path.isfile(args.genome):
            args.genome = os.path.abspath(args.genome)

    if args.emailAddress:
        # Must have at least an email server specified
        if args.smtpServer == "" or not args.smtpServer:
            sys.exit("Sorry, there is no SMTP server specified in defaults.yaml. Please specify one with --smtpServer")
        if args.emailSender == "" or not args.emailSender:
            sys.exit("Sorry, there is no email sender specified in defaults.yaml. Please specify one with --emailSender")


def commonYAMLandLogs(baseDir, workflowDir, defaults, args, callingScript):
    """
    Merge dictionaries, write YAML files, construct the snakemake command
    and create the DAG
    """
    workflowName = os.path.basename(callingScript)
    os.makedirs(args.outdir, exist_ok=True)

    if isinstance(args.snakemakeOptions, list):
        args.snakemakeOptions = ' '.join(args.snakemakeOptions)

    # save to configs.yaml in outdir
    config = defaults
    config.update(vars(args))  # This allows modifications of args after handling a user config file to still make it to the YAML given to snakemake!
    write_configfile(os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)), config)

    # merge cluster config files: 1) global one, 2) workflow specific one, 3) user provided one
    cfg = load_configfile(os.path.join(baseDir, "shared", "defaults.yaml"), False, "defaults")
    if os.path.isfile(os.path.join(baseDir, cfg['clusterConfig'])):
        cluster_config = load_configfile(os.path.join(baseDir, cfg['clusterConfig']), False)
    else:
        cluster_config = load_configfile(os.path.join(cfg['clusterConfig']), False)
    cluster_config = merge_dicts(cluster_config, load_configfile(os.path.join(workflowDir, "cluster.yaml"), False), )

    if args.clusterConfigFile:
        user_cluster_config = load_configfile(args.clusterConfigFile, False)
        cluster_config = merge_dicts(cluster_config, user_cluster_config)  # merge/override variables from user_config.yaml
    # Ensure the cluster log directory exists
    if re.search("\\{snakePipes_cluster_logDir\\}", cluster_config["snakemake_cluster_cmd"]):
        if "snakePipes_cluster_logDir" in cluster_config:
            if re.search("^/", cluster_config["snakePipes_cluster_logDir"]):
                os.makedirs(cluster_config["snakePipes_cluster_logDir"], exist_ok=True)
            else:
                os.makedirs(os.path.join(args.outdir, cluster_config["snakePipes_cluster_logDir"]), exist_ok=True)
            cluster_config["snakemake_cluster_cmd"] = re.sub("\\{snakePipes_cluster_logDir\\}", cluster_config["snakePipes_cluster_logDir"], cluster_config["snakemake_cluster_cmd"])
        else:
            sys.exit("\nPlease provide a key 'snakePipes_cluster_logDir' and value in the cluster configuration file!\n")
    write_configfile(os.path.join(args.outdir, '{}.cluster_config.yaml'.format(workflowName)), cluster_config)

    # Save the organism YAML file as {PIPELINE}_organism.yaml
    if workflowName != "preprocessing":
        orgyaml = os.path.join(baseDir, cfg['organismsDir'], "{}.yaml".format(args.genome))
        if not os.path.isfile(orgyaml):
            if os.path.isfile(os.path.join(cfg['organismsDir'], "{}.yaml".format(args.genome))):
                orgyaml = os.path.join(cfg['organismsDir'], "{}.yaml".format(args.genome))
            else:
                orgyaml = args.genome
        organismYAMLname = os.path.join(args.outdir, "{}_organism.yaml".format(workflowName))
        if workflowName != "createIndices" and os.path.abspath(organismYAMLname) != os.path.abspath(orgyaml):
            shutil.copyfile(orgyaml, organismYAMLname)

    if args.keepTemp:
        args.snakemakeOptions += " --notemp"

    snakemake_cmd = """
                    TMPDIR={tempDir}
                    UTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
                    XDG_CACHE_HOME=$UTEMP TMPDIR={tempDir} PYTHONNOUSERSITE=True snakemake {snakemakeOptions} --latency-wait {latency_wait} --snakefile {snakefile} --jobs {maxJobs} --directory {workingdir} --configfile {configFile} --keep-going --use-conda --conda-prefix {condaEnvDir}
                    """.format(latency_wait=cluster_config["snakemake_latency_wait"],
                               snakefile=os.path.join(workflowDir, "Snakefile"),
                               maxJobs=args.maxJobs,
                               workingdir=args.workingdir,
                               snakemakeOptions=str(args.snakemakeOptions or ''),
                               tempDir=cfg["tempDir"],
                               configFile=os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)),
                               condaEnvDir=cfg["condaEnvDir"]).split()

    if args.verbose:
        snakemake_cmd.append("--printshellcmds")

    if not args.local:
        snakemake_cmd += ["--cluster-config",
                          os.path.join(args.outdir, '{}.cluster_config.yaml'.format(workflowName)),
                          "--cluster", "'" + cluster_config["snakemake_cluster_cmd"], "'"]
    return " ".join(snakemake_cmd)


def print_DAG(args, snakemake_cmd, callingScript, defaults):
    if args.createDAG:
        config = defaults
        config.update(vars(args))
        workflowName = os.path.basename(callingScript)
        oldVerbose = config['verbose']
        config['verbose'] = False
        write_configfile(os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)), config)
        DAGproc = subprocess.Popen(snakemake_cmd + " --rulegraph ", stdout=subprocess.PIPE, shell=True)
        subprocess.check_call("dot -Tpdf -o{}/{}_pipeline.pdf".format(args.outdir, workflowName), stdin=DAGproc.stdout, shell=True)
        config['verbose'] = oldVerbose
        write_configfile(os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)), config)


def logAndExport(args, workflowName):
    """
    Set up logging
    """
    # Write snakemake_cmd to log file
    fnames = glob.glob(os.path.join(args.outdir, '{}_run-[0-9]*.log'.format(workflowName)))
    if len(fnames) == 0:
        n = 1  # no matching files, this is the first run
    else:
        fnames.sort(key=os.path.getctime)
        n = int(fnames[-1].split("-")[-1].split(".")[0]) + 1  # get new run number
    # append the new run number to the file name
    logfile_name = "{}_run-{}.log".format(workflowName, n)

    return logfile_name


def runAndCleanup(args, cmd, logfile_name):
    """
    Actually run snakemake. Kill its child processes on error.
    Also clean up when finished.
    """
    if args.verbose:
        print("\n{}\n".format(cmd))

    # write log file
    f = open(os.path.join(args.outdir, logfile_name), "w")
    f.write(" ".join(sys.argv) + "\n\n")
    f.write(cmd + "\n\n")

    # Run snakemake, stderr -> stdout is needed so readline() doesn't block
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if args.verbose:
        print("PID:", p.pid, "\n")

    while p.poll() is None:
        stdout = p.stdout.readline(1024)
        if stdout:
            sys.stdout.write(stdout.decode('utf-8'))
            f.write(stdout.decode('utf-8'))
            sys.stdout.flush()
            f.flush()
    # This avoids the race condition of p.poll() exiting before we get all the output
    stdout = p.stdout.read()
    if stdout:
        sys.stdout.write(stdout.decode('utf-8'))
        f.write(stdout.decode('utf-8'))
    f.close()

    # Exit with an error if snakemake encountered an error
    if p.returncode != 0:
        sys.stderr.write("Error: snakemake returned an error code of {}, so processing is incomplete!\n".format(p.returncode))
        if args.emailAddress:
            sendEmail(args, p.returncode)
        sys.exit(p.returncode)
    else:
        Path(
            os.path.join(args.outdir, "{}_snakePipes.done".format(logfile_name.split('_')[0]))
        ).touch()
        if os.path.exists(os.path.join(args.outdir, ".snakemake")):
            shutil.rmtree(os.path.join(args.outdir, ".snakemake"), ignore_errors=True)

    # Send email if desired
    if args.emailAddress:
        sendEmail(args, 0)


def predict_chip_dict(wdir, input_pattern_str, bamExt, fromBAM=None):
    """
    Predict a chip_dict from set of bam files
    ChIP input/control samples are identified from input_pattern (default: 'input')
    for each sample then the best input sample (by fuzzywuzzy score) is selected
    chip_dict is written as yaml to workflow workingdir
    predicts whether a sample is broad or narrow based on histone mark pattern
    """
    pat = "|".join(re.split(',| |\\||;', input_pattern_str))
    input_pat = r".*(" + pat + ")"
    clean_pat = r"" + pat + ""
    pat1 = re.compile(clean_pat, re.IGNORECASE)

    if fromBAM:
        infiles = sorted(glob.glob(os.path.join(fromBAM, '*' + bamExt)))
    else:
        infiles = sorted(glob.glob(os.path.join(wdir, 'filtered_bam/', '*.bam')))
    samples = get_sample_names_bam(infiles, bamExt)

    chip_dict_pred = {}
    chip_dict_pred["chip_dict"] = {}
    print("---------------------------------------------------------------------------------------")
    print("Predict Chip-seq sample configuration")
    print("---------------------------------------------------------------------------------------")
    print("\nSearch for Input/control samples...")

    input_samples = set([])
    for i in samples:
        if re.match(input_pat, i, re.IGNORECASE):
            print("...found: ", i)
            input_samples.add(i)

    print("\nTry to find corresponding ChIP samples...")

    for i in samples:
        if i in input_samples:
            continue

        print("\n sample: ", i,)
        matches_sim = {}
        for j in input_samples:
            c_clean = pat1.sub("", j)
            sim1 = fuzz.ratio(c_clean, i) + fuzz.partial_ratio(c_clean, i) + fuzz.token_sort_ratio(c_clean, i) + fuzz.token_set_ratio(c_clean, i)
            matches_sim[j] = sim1 / 4

        sim = 0
        final_matches = set([])
        for key, value in sorted(matches_sim.items(), key=lambda k: (k[1], k[0]), reverse=True):
            if value >= sim:
                final_matches.add(key)
                print("   top matching input sample by score: %s = %s" % (key, value))
                sim = value

        tmp = ':'.join(list(final_matches))

        if len(final_matches) > 1:
            tmp = "__PLEASE_SELECT_ONLY_ONE_CONTROL__:" + tmp
        elif len(final_matches) == 0:
            print("No control sample found!")

        chip_dict_pred["chip_dict"][i] = {}
        chip_dict_pred["chip_dict"][i]['control'] = tmp
        if re.match(".*(H3K4me1|H3K36me3|H3K9me3|H3K27me3).*", i, re.IGNORECASE):
            chip_dict_pred["chip_dict"][i]['broad'] = True
        else:
            chip_dict_pred["chip_dict"][i]['broad'] = False

    outfile = os.path.join(wdir, "chip_seq_sample_config.PREDICTED.yaml")
    write_configfile(outfile, chip_dict_pred)
    print("---------------------------------------------------------------------------------------")
    print("Chip-seq sample configuration is written to file ", outfile)
    print("Please check and modify this file - this is just a guess! Then run the workflow with it.")
    print("---------------------------------------------------------------------------------------")


def writeTools(usedEnvs, wdir, workflowName, maindir):
    outfile = os.path.join(wdir, workflowName + "_tools.txt")
    with open(outfile, 'w') as f:
        for item in usedEnvs:
            dependencies = False
            for line in open(os.path.join(maindir, "shared", "rules", item), 'r'):
                if line.split(":")[0] == "dependencies":
                    dependencies = True
                elif dependencies is True:
                    f.write(line)
    f.close()


def copySampleSheet(sampleSheet, wdir):
    if os.path.isfile(sampleSheet) and os.path.exists(wdir):
        bname = os.path.basename(sampleSheet)
        try:
            shutil.copyfile(sampleSheet, os.path.join(wdir, bname))
        except Exception as err:
            print("Unexpected error:\n{}".format(err))
            raise
