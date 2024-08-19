#!/usr/bin/env python3

__description__ = """
MPI-IE workflow for preprocessing

usage example:
    preprocessing -i input-dir -o output-dir --optDedupDist 2500 --clumpifyOptions
"""


import argparse
import os
import sys
import glob
import textwrap
import snakePipes.common_functions as cf
import snakePipes.parserCommon as parserCommon


def parse_args(defaults={"verbose": False, "configFile": None,
                         "clusterConfigFile": None, "maxJobs": 5,
                         "snakemakeOptions": "--use-conda", "tempDir": None,
                         "fastqc": False, "optDedupDist": 0,
                         "clumpifyOptions": "dupesubs=0 qin=33 markduplicates=t optical=t",
                         "clumpifyMemory": "30G", "sampleSheet": None,
                         "ext": ".fastq.gz", "reads": ["_R1", "_R2"],
                         "UMIDedupSep": "_", "UMIBarcode": False,
                         "bcPattern": ""}):
    """
    Parse arguments from the command line.
    """
    mainArgs = parserCommon.mainArguments(defaults, preprocessing=True)

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        parents=[mainArgs],
        add_help=False
    )

    # Workflow options
    optional = parser.add_argument_group('Options')
    parserCommon.commonOptions(optional, defaults, preprocessing=True)

    optional.add_argument("--optDedupDist",
                          type=int,
                          help="The maximum distance between clusters to mark one as "
                          "an optical duplicate of the other. Setting this to 0 will "
                          "disable optical deduplication, which is only needed on "
                          "patterned flow cells (NextSeq, NovaSeq, HiSeq 3000/4000, "
                          "etc.). Common values are: 2500 (HiSeq 3000/4000), 40 "
                          "(NextSeq) or 12000 (NovaSeq). (default: '%(default)s')",
                          default=defaults["optDedupDist"])

    optional.add_argument("--clumpifyOptions",
                          help="Options passed to clumpify, which should generally "
                          "NOT be changed. The exception to this is with NextSeq runs, "
                          "where 'spany=t adjacent=t' should be ADDED. (default: '%(default)s)",
                          default=defaults["clumpifyOptions"])

    optional.add_argument("--clumpifyMemory",
                          help="This controls how much memory clumpify is instructed "
                          "to use, in GB. This may need to be increased if "
                          "samples are particularly large or there are MANY optical "
                          "duplicates. This is passed to clumpify (e.g., as "
                          "-Xmx30G). You may additionally need to instruct your "
                          "scheduler about the per-core memory usage (e.g., in cluster.yaml). "
                          "(default: '%(default)s')",
                          default=defaults["clumpifyMemory"])

    mutex = optional.add_mutually_exclusive_group()
    mutex.add_argument("--sampleSheet",
                       help="Information on samples (required for merging across lanes); see "
                            "'docs/content/preprocessing_sampleSheet.example.tsv' for an example."
                            " The first set of columns should hold the current sample names"
                            " (excluding things like _R1.fastq.gz or _1.fastq.gz) while the"
                            " second holds the name that the final sample should have (again,"
                            " excluding things like _R1.fastq.gz or _1.fastq.gz)."
                            " (default: '%(default)s')",
                       default=defaults["sampleSheet"])

    mutex.add_argument("--automaticIlluminaMerging",
                       help="Instead of manually specifying a sample sheet that can be used to "
                            "merge samples, this option assumes that fastq files in the input "
                            "directory follow the default output of bcl2fastq from Illumina and "
                            "should be merged across lanes. An example of this would be that "
                            "files named Sample1_S1_L001_R1_001.fastq.gz and "
                            "Sample1_S1_L002_R1_001.fastq.gz would be merged to "
                            "Sample1_R1.fastq.gz. This option CAN NOT be combined with "
                            "--sampleSheet. It is assumed that standard R1/R2 read designators "
                            "are being used and that all files end in .fastq.gz",
                       action='store_true')

    return parser


def generateSampleSheet(indir, reads, ext, outdir):
    ofile = os.path.join(outdir, "generatedSampleSheet.txt")
    of = open(ofile, "w")
    for fname in sorted(glob.glob(os.path.join(indir, "*_S[0-9]*_L[0-9]*_R?_001.fastq.gz"))):
        fname = os.path.basename(fname)
        sampleName = fname[:-9].split("_")
        readDesignator = sampleName[-2]
        sampleName = "_".join(sampleName[:-4])
        if readDesignator == "R1":
            of.write("{}\t{}\t{}\n".format(fname, reads[0], sampleName))
        else:
            of.write("{}\t{}\t{}\n".format(fname, reads[1], sampleName))
    of.close()
    return ofile


def main():
    baseDir, workflowDir, defaults = cf.setDefaults(os.path.basename(__file__))

    # get command line arguments
    parser = parse_args(defaults)
    args = parser.parse_args()
    args, defaults = cf.handleUserArgs(args, defaults, parse_args)

    # we also add these paths to config, although we don't use them in the Snakefile
    args.baseDir = baseDir

    # Common arguments
    cf.checkCommonArguments(args, baseDir, outDir=True, preprocessing=True)

    # Generate a sample sheet if needed
    if args.automaticIlluminaMerging:
        args.sampleSheet = generateSampleSheet(args.indir, args.reads, args.ext, args.outdir)

    # Handle YAML and log files
    snakemake_cmd = cf.commonYAMLandLogs(baseDir, workflowDir, defaults, args, __file__)
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # Run everything
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

    #CreateDAG
    cf.plot_DAG(args,snakemake_cmd, __file__,defaults)


if __name__ == "__main__":
    main()
