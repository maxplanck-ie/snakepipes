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


def set_env_yamls():
    """
    This defines the global variables describing where the conda env yaml files are
    """
    return {'CONDA_SHARED_ENV': 'envs/shared.yaml',
            'CONDA_CREATE_INDEX_ENV': 'envs/createIndices.yaml',
            'CONDA_RNASEQ_ENV': 'envs/rna_seq.yaml',
            'CONDA_scRNASEQ_ENV': 'envs/sc_rna_seq.yaml',
            'CONDA_DNA_MAPPING_ENV': 'envs/dna_mapping.yaml',
            'CONDA_CHIPSEQ_ENV': 'envs/chip_seq.yaml',
            'CONDA_ATAC_ENV': 'envs/atac_seq.yaml',
            'CONDA_HIC_ENV': 'envs/hic.yaml',
            'CONDA_WGBS_ENV': 'envs/wgbs.yaml',
            'CONDA_PY27_ENV': 'envs/python27.yaml',
            'CONDA_RMD_ENV': 'envs/rmarkdown.yaml'}


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


def load_configfile(configfile, verbose, info='Config'):
    with open(configfile, "r") as f:
        config = yaml.load(f)

    config = sanity_dict_clean(config)

    if verbose:
        print("\n--- " + info + " ---------------------------------------------------------------------")
        print("config file: {}".format(configfile))
        for k, v in sorted(config.items()):
            print("{}: {}".format(k, v))
        print("-" * 80, "\n")
    return config


def write_configfile(configfile, config):
    with open(configfile, 'w') as f:
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


def load_organism_data(genome, maindir, verbose):
    if os.path.isfile(os.path.join(maindir, "shared", "organisms", genome + ".yaml")):
        organism = load_configfile(os.path.join(maindir, "shared", "organisms", genome + ".yaml"), verbose, "Genome")
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
        elif x.endswith(reads[1]):
            x = x[:-l1]
        else:
            continue
        s.add(x)
    return sorted(list(s))


def get_sample_names_bam(infiles, bam_ext):
    """
    Get sample names without file extensions
    """
    s = []
    for x in infiles:
        x = os.path.basename(x).replace(bam_ext, "")
        s.append(x)
    return sorted(list(set(s)))


def is_paired(infiles, ext, reads):
    """
    Check for paired-end input files
    """
    paired = False
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
    if infiles_dic and max([len(x) for x in infiles_dic.values()]) == 2:
        paired = True
    # TODO: raise exception if single-end and paired-end files are mixed
    return paired


def get_fragment_length(infile, sampleName):
    """
    Return median insert size from a metrics file created by
    deeptools bamPEFragmentSize.
    Read the 37 column text file, grep the row corresponding to the sample and
    return the entry from 6th column (Frag. Len. Median)
    """
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("filtered_bam/{}".format(sampleName)):
                try:
                    median = line.split()[5]
                    return int(float(median))
                except TypeError:
                    print("ERROR: File", infile, "is NOT an output from bamPEFragmentSize.\n")
                    exit(1)
            else:
                pass
    # no match in infile
    print("ERROR: File", infile, "is NOT an output from bamPEFragmentSize.\n")
    exit(1)


def make_temp_dir(tempdir, fallback_dir, verbose=False):
    try:
        output = subprocess.check_output("mktemp -d -p " + tempdir + "/ tmp.snakemake.XXXXXXXX", shell=True, stderr=subprocess.STDOUT)
        temp_path = output.decode().rstrip() + "/"
    except subprocess.CalledProcessError:
        try:
            print("\nFailed to create temp dir under temp path prefix (" + tempdir + ")! Try fallback: " + fallback_dir + "/ ...")
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
                if args.strains is '':
                    print("\nError! Please specify strain ID to extract from given VCF file for Allele-specific mapping! ({})\n".format(args.VCFfile))
                    exit(1)
                else:
                    allele_mode = 'create_and_map'
            else:
                print("\nError! Please specify either VCF file or SNP file for Allele-specific mapping! \n")
                exit(1)
        # If SNP file is present, check whether genome index also exists
        elif not os.path.exists(os.path.dirname(args.Nmasked_index)):
            print("\nError! Please specify an n-masked index file for Allele-specific mapping! \n")
            exit(1)
        else:
            allele_mode = 'map_only'
    else:
        allele_mode = None
    return allele_mode


def cleanLogs(d):
    """
    Remove all empty log files, both in cluster_logs/ and */logs/
    """
    for f in glob.glob(os.path.join(d, "cluster_logs", "*")):
        s = os.stat(f)
        if s.st_size == 0:
            os.remove(f)
    for f in glob.glob(os.path.join(d, "*", "logs", "*")):
        s = os.stat(f)
        if s.st_size == 0:
            os.remove(f)


def check_sample_info_header(sample_info_file):
    """
    return True in case sample info file contains column names 'name' and 'condition'
    """
    ret = subprocess.check_output("cat " + sample_info_file + " | head -n1",
                                  shell=True).decode()

    if "name" in ret.split() and "condition" in ret.split():
        return True
    else:
        return False


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


def sendEmail(args, returnCode):
    """
    Try to send an email to the user. Errors must be non-fatal.
    """
    try:
        import smtplib
        from email.message import EmailMessage
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


def checkCommonArguments(args, baseDir, outDir=False, createIndices=False):
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
            if os.path.exists(args.workingdir):
                args.workingdir = os.path.abspath(args.workingdir)
            else:
                sys.exit("\nError! Working-dir (-d) dir not found! ({})\n".format(args.workingdir))
            args.outdir = args.workingdir
    args.cluster_logs_dir = os.path.join(args.outdir, "cluster_logs")
    # 2. Config file
    if args.configfile and not os.path.exists(args.configfile):
        sys.exit("\nError! Provided configfile (-c) not found! ({})\n".format(args.configfile))
    # 3. Sample info file
    if 'sample_info' in args and args.sample_info:
        if os.path.exists(os.path.abspath(args.sample_info)):
            args.sample_info = os.path.abspath(args.sample_info)
        else:
            sys.exit("\nSample info file not found! (--DB {})\n".format(args.sample_info))
        if not check_sample_info_header(args.sample_info):
            sys.exit("ERROR: Please use 'name' and 'condition' as column headers in sample info file! ({})\n".format(args.sample_info))
    # 4. get abspath from user provided genome/organism file
    if not createIndices:
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
    snakemake_path = os.path.dirname(os.path.abspath(callingScript))

    # merge configuration dicts
    config = defaults   # 1) form defaults.yaml
    if args.configfile:
        user_config = load_configfile(args.configfile, False)
        config = merge_dicts(config, user_config)  # 2) from user_config.yaml
    config_wrap = config_diff(vars(args), defaults)  # 3) from wrapper parameters
    config = merge_dicts(config, config_wrap)

    # Ensure the log directory exists
    os.makedirs(args.cluster_logs_dir, exist_ok=True)

    # save to configs.yaml in outdir
    write_configfile(os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)), config)

    # merge cluster config files: 1) global one, 2) workflow specific one, 3) user provided one
    cluster_config = load_configfile(os.path.join(baseDir, "shared/cluster.yaml"), False)
    cluster_config = merge_dicts(cluster_config, load_configfile(os.path.join(workflowDir, "cluster.yaml"), False), )

    if args.cluster_configfile:
        user_cluster_config = load_configfile(args.cluster_configfile, False)
        cluster_config = merge_dicts(cluster_config, user_cluster_config)  # merge/override variables from user_config.yaml
    write_configfile(os.path.join(args.outdir, '{}.cluster_config.yaml'.format(workflowName)), cluster_config)

    # Save the organism YAML file as {PIPELINE}_organism.yaml
    orgyaml = os.path.join(baseDir, "shared/organisms/{}.yaml".format(args.genome))
    if not os.path.isfile(orgyaml):
        orgyaml = args.genome
    organismYAMLname = os.path.join(args.outdir, "{}_organism.yaml".format(workflowName))
    if os.path.abspath(organismYAMLname) != os.path.abspath(orgyaml):
        shutil.copyfile(orgyaml, organismYAMLname)

    if args.notemp:
        args.snakemake_options += " --notemp"

    snakemake_cmd = """
                    {snakemake} {snakemake_options} --latency-wait {latency_wait} --snakefile {snakefile} --jobs {max_jobs} --directory {workingdir} --configfile {configfile} --keep-going
                    """.format(snakemake=os.path.join(snakemake_path, "snakemake"),
                               latency_wait=cluster_config["snakemake_latency_wait"],
                               snakefile=os.path.join(workflowDir, "Snakefile"),
                               max_jobs=args.max_jobs,
                               workingdir=args.workingdir,
                               snakemake_options=str(args.snakemake_options or ''),
                               configfile=os.path.join(args.outdir, '{}.config.yaml'.format(workflowName))).split()

    # Produce the DAG if desired
    if args.createDAG:
        oldVerbose = config['verbose']
        config['verbose'] = False
        write_configfile(os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)), config)
        DAGproc = subprocess.Popen(snakemake_cmd + ['--rulegraph'], stdout=subprocess.PIPE)
        _ = open("{}/{}_pipeline.pdf".format(args.outdir, workflowName), "wb")
        subprocess.check_call(["dot", "-Tpdf"], stdin=DAGproc.stdout, stdout=_)
        _.close()
        config['verbose'] = oldVerbose
        write_configfile(os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)), config)

    if args.verbose:
        snakemake_cmd.append("--printshellcmds")

    if not args.local:
        snakemake_cmd += ["--cluster-config",
                          os.path.join(args.outdir, '{}.cluster_config.yaml'.format(workflowName)),
                          "--cluster", "'" + cluster_config["snakemake_cluster_cmd"],
                          args.cluster_logs_dir, "--name {rule}.snakemake'"]

    return snakemake_cmd


def logAndExport(args, workflowName):
    """
    Set up logging and exports (TMPDIR)
    """
    # Write snakemake_cmd to log file
    fnames = glob.glob(os.path.join(args.outdir, '{}_run-[0-9*].log'.format(workflowName)))
    if len(fnames) == 0:
        n = 1  # no matching files, this is the first run
    else:
        fnames.sort(key=os.path.getctime)
        n = int(fnames[-1].split("-")[-1].split(".")[0]) + 1  # get new run number
    # append the new run number to the file name
    logfile_name = "{}_run-{}.log".format(workflowName, n)

    # create local temp dir and add this path to environment as $TMPDIR variable
    # on SLURM: $TMPDIR is set, created and removed by SlurmEasy on cluster node
    temp_path = make_temp_dir(args.tempdir, args.outdir)
    snakemake_exports = ["export", "TMPDIR='{}'".format(temp_path), "&&"]

    return snakemake_exports, logfile_name, temp_path


def runAndCleanup(args, cmd, logfile_name, temp_path):
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

    # remove temp dir
    if (temp_path != "" and os.path.exists(temp_path)):
        shutil.rmtree(temp_path, ignore_errors=True)
        if args.verbose:
            print("Temp directory removed ({})!\n".format(temp_path))

    # Send email if desired
    if args.emailAddress:
        sendEmail(args, 0)


def predict_chip_dict(wdir):
    """
    Predict a chip_dict from bam files under filtered_bam/ from DNA-mapping workflow
    ChIP input/control samples are identified from pattern 'input' (case ignored)
    chip_dict is written as yaml to current workflow workingdir
    predicts whether a sample is broad or narrow based on histone mark pattern
    """
    pat1 = re.compile(r"input.*$", re.IGNORECASE)
    pat2 = re.compile(r"^.*input", re.IGNORECASE)
    infiles = sorted(glob.glob(os.path.join(wdir, 'filtered_bam/', '*.bam')))
    samples = get_sample_names(infiles, ".filtered.bam", ['', ''])

    chip_dict_pred = {}
    chip_dict_pred["chip_dict"] = {}
    print("---------------------------------------------------------------------------------------")
    print("Predict Chip-seq sample configuration")
    print("---------------------------------------------------------------------------------------")
    print("\nSearch for Input/control samples...")

    input_samples = set([])
    for i in samples:
        if re.match(r".*input.*", i, re.IGNORECASE):
            print("...found: ", i)
            input_samples.add(i)

    print("\nTry to find corresponding ChIP samples...")
    final_matches = set()
    for i in samples:
        if i in input_samples:
            continue

        print("\n sample: ", i)

        prefix_matches = set([])
        suffix_matches = set([])

        for j in input_samples:
            c_prefix = pat1.sub("", j)
            c_suffix = pat2.sub("", j)

            if re.match(r"^" + c_prefix + ".*", i, re.IGNORECASE):
                prefix_matches.add(j)
            if re.match(r".*" + c_suffix + "$", i, re.IGNORECASE):
                suffix_matches.add(j)

        final_matches = set([])

        if len(prefix_matches) > 0:
            final_matches = prefix_matches

        if len(suffix_matches) > 0 and (len(prefix_matches) == 0 or len(suffix_matches) < len(prefix_matches)):
            final_matches = suffix_matches

        if len(prefix_matches) == len(suffix_matches) and len(prefix_matches) > 0:
            final_matches = set(prefix_matches).update(suffix_matches)

        tmp = ':'.join(list(final_matches))
        print("   pref:", prefix_matches, " suf:", suffix_matches, " final:", tmp)

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

    write_configfile(os.path.join(wdir, "chip_seq_sample_config.yaml"), chip_dict_pred)
    print("---------------------------------------------------------------------------------------")
    print("Chip-seq sample configuration is written to file ", os.path.join(wdir, "chip_seq_sample_config.yaml"))
    print("Please check and modify this file - this is just a guess! Then run the workflow with it.")
    print("---------------------------------------------------------------------------------------")
    sys.exit(0)
