#!/usr/bin/env python
import sys
import argparse
import subprocess
import snakePipes
import os
import yaml
import glob
import hashlib
import shutil
import snakePipes.common_functions as cof
from importlib.metadata import version


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Setup and information script for snakePipes",
        usage="$ snakePipes info",
    )

    subparsers = parser.add_subparsers(title="Commands", dest="command")

    subparsers.add_parser(
        "info", help="Print the location of the various yaml files"
    )

    createEnvsParser = subparsers.add_parser(
        "createEnvs",
        help="Create or update conda enviroments according to the "
        "workflow-specific yaml files. Note that changing the snakemakeOptions: "
        "option will result in ALL conda environments being recreated.",
    )

    subparsers.add_parser(
        "envInfo",
        help="Prints the location in which each conda environment is actually stored.",
    )


    createEnvsParser.add_argument(
        "--autodetectCondaEnvDir",
        action="store_true",
        help="If specified, this will set condaEnvDir to system conda prefix,"
       "and will overwrite the condaEnvDir entry in defaults.yaml ."
       "Use with caution."
    )


    createEnvsParser.add_argument(
        "--only",
        nargs="+",
        help="If specified, a space-separated list of environments to create. "
        "This should typically only be done for testing purposes. The "
        "possible environments are: {}".format(cof.set_env_yamls().keys()),
    )

    createEnvsParser.add_argument(
        "--info",
        "-i",
        action="store_true",
        help="Only print the environments that would be created, don't actually create them.",
    )

    createEnvsParser.add_argument(
        "--noSitePackages",
        action="store_true",
        help="Prevent conda from looking at anything installed in a users home "
        "directory. While this is convenient to ensure that ONLY packages "
        "installed by snakePipes are used, it means that you will "
        "occasionally get scary-looking warning messages when you try to "
        "create new conda environments.",
    )

    subparsers.add_parser(
        "flushOrganisms",
        help="Flush all installed organism YAML files. This is only advisable when performing a clean installation.",
    )

    subparsers.add_parser("version", help="Print the snakePipes version")

    baseDir = os.path.dirname(snakePipes.__file__)
    defaults = cof.load_configfile(
        os.path.join(baseDir, "shared", "defaults.yaml"), False, "defaults"
    )
    configParser = subparsers.add_parser(
        "config", help="Update snakePipes-wide values in defaults.yaml"
    )

    configParser.add_argument(
        "--configMode",
        help="Choose between manual argument setup and recycling a config file from previous installation (Default: %(default)s)",
        choices=["manual", "recycle"],
        default=defaults["configMode"],
    )

    configParser.add_argument(
        "--oldConfig",
        help="Provide an existing config file with the setup used in the previous installation. (Default: %(default)s)",
        default=defaults["oldConfig"],
    )

    configParser.add_argument(
        "--snakemakeOptions",
        help="Update the options given to snakeMake. You MUST include --use-conda "
        "and an appropriate --conda-prefix if you change this! (Default: %(default)s)",
        default=defaults["snakemakeOptions"],
    )

    configParser.add_argument(
        "--condaEnvDir",
        help="If specified, use this as the base directory for the "
        "created conda environments. This will ignore what is already "
        "in the workflow-specific yaml files and where conda is installed.",
    )

    configParser.add_argument(
        "--organismsDir",
        help="The directory where global organism YAML files are to be stored. Both "
        "absolute and relative paths are supported. In the latter case the "
        "path is then relative to the snakePipes installation directory. (Default: %(default)s)",
        default=defaults["organismsDir"],
    )

    configParser.add_argument(
        "--clusterConfig",
        help="The YAML file containing the snakeMake cluster command and global "
        "memory settings. Both absolute and relative paths are supported. "
        "In the latter case the path is then relative to the snakePipes "
        "installation directory. (Default: %(default)s)",
        default=defaults["clusterConfig"],
    )

    configParser.add_argument(
        "--tempDir",
        help="A custom directory where temporary files should be written. This "
        "is ideally locally attached to your cluster nodes. "
        "(Default: %(default)s)",
        default=defaults["tempDir"],
    )

    configParser.add_argument(
        "--noToolsVersion",
        dest="toolsVersion",
        help="By default, tool versions are printed to a workflow-specific file. Specifying this disables that behavior.",
        action="store_false",
    )

    email = configParser.add_argument_group(
        title="Email/SMTP options",
        description="These options are only used if/when --emailAddress is used.",
    )
    email.add_argument(
        "--smtpServer",
        help="SMTP server address. (Default: %(default)s)",
        default=defaults["smtpServer"],
    )

    email.add_argument(
        "--smtpPort",
        type=int,
        help="The port on the SMTP server to use. A value of 0 will use the default SMTP port. (Default: %(default)s)",
        default=defaults["smtpPort"],
    )

    email.add_argument(
        "--onlySSL",
        action="store_true",
        help="If specified, only use SSL-enabled connections.",
    )

    email.add_argument(
        "--emailSender",
        help="The email address used to send emails. (Default: %(default)s)",
        default=defaults["emailSender"],
    )

    email.add_argument(
        "--smtpUsername",
        help="For SMTP servers requiring a login, the username to use. (Default: %(default)s)",
        default=defaults["smtpUsername"],
    )

    email.add_argument(
        "--smtpPassword",
        help="For SMTP servers requiring a login, the password to use. Note that this is stored in clear text! (Default: %(default)s)",
        default=defaults["smtpPassword"],
    )

    return parser


def flushOrganisms():
    """
    Remove all organism YAML files.
    """
    baseDir = os.path.dirname(snakePipes.__file__)
    for f in glob.glob(os.path.join(baseDir, "shared/organisms/*.yaml")):
        os.remove(f)


def info():
    """
    Print the locations of EVERY yaml file. Break these up a bit so it's clear what they actually belong to. Print path to tempDir and check that it exists.
    """
    baseDir = os.path.dirname(snakePipes.__file__)
    cfg = cof.load_configfile(
        os.path.join(baseDir, "shared", "defaults.yaml"), False, "defaults"
    )

    # Organism yaml files
    print("Organism YAML files:")
    orgDir = cfg["organismsDir"]
    if not os.path.exists(orgDir):
        orgDir = os.path.join(baseDir, orgDir)
    for f in glob.glob(os.path.join(orgDir, "*.yaml")):
        print("    {}".format(f))

    # defaults.yaml under shared
    print(
        "\nThe workflow-generic defaults.yaml file is:\n    {}".format(
            os.path.join(baseDir, "shared/defaults.yaml")
        )
    )

    # cluster.yaml
    clusterConfig = cfg["clusterConfig"]
    if not os.path.isfile(clusterConfig):
        clusterConfig = os.path.join(baseDir, clusterConfig)
    print(
        "\nThe default cluster.yaml file. Its defaults are overridden by the per-workflow cluster.yaml files:\n    {}".format(
            clusterConfig
        )
    )
    print("\nWorkflow-specific cluster.yaml and defaults.yaml files are in:")
    for f in glob.glob(os.path.join(baseDir, "workflows/*/cluster.yaml")):
        print("    {}".format(os.path.dirname(f)))

    # tempDir
    tempDir = cfg["tempDir"]
    msg = ["\nTemp dir under {} ".format(tempDir)]
    if os.path.isdir(tempDir):
        msg.append("exists and will be used.")
    else:
        msg.append("does not exist and /tmp will be used instead.")
    print("".join(msg))


def envInfo():
    """
    For each environment yaml file print where its conda env is actually located
    """
    baseDir = os.path.dirname(snakePipes.__file__)

    f = open(os.path.join(baseDir, "shared/defaults.yaml"))
    cf = yaml.load(f, Loader=yaml.FullLoader)
    f.close()

    condaEnvDir=cf["condaEnvDir"]


    for env in cof.set_env_yamls().values():
        # Hash the file ala snakemake
        md5hash = hashlib.md5()
        md5hash.update(condaEnvDir.encode())
        f = open(os.path.join(baseDir, "shared/rules", env), "rb")
        md5hash.update(f.read())
        f.close()
        h = md5hash.hexdigest()

        print("{} is in:\n    {}\n".format(env, os.path.join(condaEnvDir, h)))


def fixSitePy(envPath):
    """
    We would really like to prevent any snakePipes environment from using the user site packages.
    """
    for fname in glob.glob("{}/lib/python*/site.py".format(envPath)):
        f = open(fname).read()
        lines = f.split("\n")
        lines = [
            line
            if not line.startswith("ENABLE_USER_SITE")
            else "ENABLE_USER_SITE = False"
            for line in lines
        ]
        f = open(fname, "w")
        f.write("\n".join(lines))
        f.close()

        cmd = [os.path.join(envPath, "bin", "python"), "-m", "compileall", fname]
        subprocess.check_call(cmd)


def createCondaEnvs(args):
    """
    Create all of the conda environments
    """
    baseDir = os.path.dirname(snakePipes.__file__)

    f = open(os.path.join(baseDir, "shared/defaults.yaml"))
    cf = yaml.load(f, Loader=yaml.FullLoader)
    f.close()
    condaEnvDir=cf["condaEnvDir"]
    condaDirUse=condaEnvDir

    if args.autodetectCondaEnvDir:
        condaDirUse=detectCondaDir()
    # rewrite defaults.yaml
        cof.write_configfile(os.path.join(baseDir, "shared/defaults.yaml"), cf)

    for envName, env in cof.set_env_yamls().items():
        if args.only is not None and envName not in args.only:
            continue
        # Hash the file ala snakemake
        md5hash = hashlib.md5()
        md5hash.update(condaDirUse.encode())
        f = open(os.path.join(baseDir, "shared/rules", env), "rb")
        md5hash.update(f.read())
        f.close()
        h = md5hash.hexdigest()

        sys.stderr.write(
            "Creating environment from {} in {}\n".format(
                os.path.join(baseDir, "shared/rules", env), condaDirUse
            )
        )
        cmd = [
            "mamba",
            "env",
            "create",
            "--file",
            os.path.join(baseDir, "shared/rules", env),
        ]
        cmd += ["--prefix", os.path.join(condaDirUse, h)]

        # Don't actually create the env if either --info is set or it already exists and --force is NOT set
        if not args.info:
            if not os.path.exists(os.path.join(condaDirUse, h)):
                try:
                    os.makedirs(os.path.join(condaDirUse, h), exist_ok=True)
                    subprocess.check_call(cmd)
                except:
                    # Ensure an environment is fully removed on error
                    shutil.rmtree(os.path.join(condaDirUse, h), ignore_errors=False)
                    sys.exit("There was an error when creating the environments!\n")

        # Ignore site-packages
        if args.noSitePackages and not args.info:
            fixSitePy(os.path.join(condaDirUse, h))

    # Ignore site-packages in this env
    if args.noSitePackages and not args.info:
        fixSitePy(rootDir)

def detectCondaDir():
    "Detect the default conda folder."
    condaDir = os.environ.get("CONDA_PREFIX")
    if "envs" in condaDir:
        condaDir = os.path.dirname(condaDir)
    else:
        condaDir = os.path.join(condaDir, "envs")
    return(condaDir)


def updateConfig(args):
    """Update the global defaults"""
    baseDir = os.path.dirname(snakePipes.__file__)
    # Load, update and rewrite the default dictionary
    currentDict = cof.load_configfile(
        os.path.join(baseDir, "shared", "defaults.yaml"), False, "Default Config"
    )

    if args.configMode == "manual":
        d = {
            "snakemakeOptions": args.snakemakeOptions,
            "condaEnvDir": args.condaEnvDir,
            "organismsDir": args.organismsDir,
            "clusterConfig": args.clusterConfig,
            "tempDir": args.tempDir,
            "smtpServer": args.smtpServer,
            "smtpPort": args.smtpPort,
            "onlySSL": args.onlySSL,
            "emailSender": args.emailSender,
            "smtpUsername": args.smtpUsername,
            "smtpPassword": args.smtpPassword,
            "toolsVersion": args.toolsVersion,
            "oldConfig": None,
            "configMode": "manual",
        }
    elif args.configMode == "recycle":
        oldConfig = args.oldConfig
        if os.path.isfile(oldConfig):
            d = cof.load_configfile(oldConfig, False, "Old Config")
            if args.organismsDir:
                od = {"organismsDir": args.organismsDir}
                d.update(od)
            if args.clusterConfig:
                od = {"clusterConfig": args.clusterConfig}
                d.update(od)
            if not currentDict.keys() & d.keys():
                sys.exit("The old and the new config have no matching keys!!!\n")
        else:
            sys.exit("Config file not found\n")
    updatedDict = cof.merge_dicts(currentDict, d)
    cof.write_configfile(os.path.join(baseDir, "shared", "defaults.yaml"), updatedDict)
    cof.load_configfile(
        os.path.join(baseDir, "shared", "defaults.yaml"), True, "Final Updated Config"
    )

def main():
    if len(sys.argv) == 1:
        sys.argv.append("--help")
    args = parse_arguments().parse_args(sys.argv[1:])
    if args.command == "info":
        info()
    elif args.command == "envInfo":
        envInfo()
    elif args.command == "flushOrganisms":
        flushOrganisms()
    elif args.command == "config":
        updateConfig(args)
    elif args.command == "version":
        _v = version("snakePipes")
        print(f"snakePipes version {_v}")
    else:
        createCondaEnvs(args)
