#!/usr/bin/env python
# functions shared across workflows ##########################################
################################################################################
import subprocess
import os
import re
import yaml


def get_snakepipes_path():

    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def convert_library_type(R_path, paired, from_library_type, from_prg, to_prg, rscript, tsv):
    """ Converts the library to e.g. from 2 (featureCounts) to RF (HISAT2) """
    if paired:
        lib_str = "PE"
    else:
        lib_str = "SE"

    cmd = ("{}Rscript {} {} {} {} {} {}".format(R_path, rscript, tsv, lib_str, from_library_type, from_prg, to_prg))
    # print("\n"+cmd)

    return subprocess.check_output(cmd, shell=True).decode()


def merge_dicts(x, y):
    z = {}
    z = x.copy()
    if y:
        z.update(y)
    return z


# this is a pure sanity fucntion to avoid obvious mailfunction during snakefile execution
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


def load_paths(pathfile, maindir, verbose):
    paths = load_configfile(pathfile, False)

    # add path to tools dir
    paths["workflow_tools"] = os.path.join(maindir, "shared", "tools")

    if verbose:
        print("\n--- paths ---------------------------------------------------------------------")
        for k, v in sorted(paths.items()):
            print("{}: {}".format(k, v))
        print("-" * 80, "\n")

    return paths


def get_sample_names(infiles, ext, reads):
    """
    Get sample names without file extensions
    """
    s = []
    for x in infiles:
        x = os.path.basename(x).replace(ext, "")
        try:
            x = x.replace(reads[0], "").replace(reads[1], "")
        except IndexError:
            pass
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
            # print(m.group())
            bname = m.group(1)
            # print(bname)
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
                    return float(median)
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
