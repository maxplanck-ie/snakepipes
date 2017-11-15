#!/usr/bin/env python3


# paired-end ONLY!


# Usage: python sample_qc_report.py sample_name.alignment_summary_metrics.txt sample_name.mark_duplicates_metrics.txt sample_name.filtered.BAM_peaks.xls sample_name.filtered.BAM_peaks.qc.txt
# python ~/git/snakemake_workflows/shared/tools/sample_qc_report_PE.py Picard_qc/AlignmentSummaryMetrics/Cor29_79me2.alignment_summary_metrics.txt Picard_qc/MarkDuplicates/Cor29_79me2.mark_duplicates_metrics.txt Picard_qc/InsertSizeMetrics/Cor29_79me2.insert_size_metrics.txt MACS2/Cor29_79me2.filtered.BAM_peaks.xls MACS2/Cor29_79me2.filtered.BAM_peaks.qc.txt


import os
import sys


# Input files ################################################################
infile_AlignmentSummaryMetrics = sys.argv[1]    # sample_name.alignment_summary_metrics.txt
infile_MarkDuplicates = sys.argv[2]             # sample_name.mark_duplicates_metrics.txt
input_InsertSizeMetrics = sys.argv[3]
try:
    infile_MACS2_xls = sys.argv[4]              # sample_name.filtered.BAM_peaks.xls (optional)
    infile_MACS2_qc_txt = sys.argv[5]           # ample_name.filtered.BAM_peaks.qc.txt (optional)
except IndexError:
    pass


# read Picard AlignmentSummaryMetrics results ################################
try:
    with open(infile_AlignmentSummaryMetrics) as f:
        lines = f.readlines()

    columns = list(filter(lambda x: x.startswith("PAIR"), lines))[0].strip().split("\t")

    # PF_READS: The number of PF reads where PF is defined as passing Illumina's filter.
    total_reads = int(columns[2])

    # PF_READS_ALIGNED: The number of PF reads that were aligned to the reference sequence.
    mapped_reads = int(columns[5])

    # PCT_PF_READS_ALIGNED: The percentage of PF reads that aligned to the reference sequence. PF_READS_ALIGNED / PF_READS
    fmapped_reads = float(columns[6])

    # fraction of high-quality mapped reads (MAPQ >=20)
    fhq_mapped_reads = 1.0 * int(columns[8]) / total_reads

    # specific for paired-end
    read_pairs = total_reads / 2
    mapped_reads_in_pairs = int(columns[16])  # READS_ALIGNED_IN_PAIRS
    mapped_pairs = mapped_reads_in_pairs / 2
    fmapped_pairs = 1.0 * mapped_pairs / read_pairs
    fmapped_singletons = 1.0 * (mapped_reads - mapped_reads_in_pairs) / total_reads

except OSError:
    exit("ERROR! Unable to read: {}".format(infile_AlignmentSummaryMetrics))


# read Picard MarkDuplicates results #########################################
try:
    with open(infile_MarkDuplicates) as f:
        lines = f.readlines()
    columns = list(filter(lambda x: x.startswith("Unknown Library"), lines))[0].strip().split("\t")

    dup_mapped_pairs = int(columns[5])
    fdup_mapped_pairs = 1.0 * dup_mapped_pairs / mapped_pairs
    dupfree_mapped_pairs = mapped_pairs - dup_mapped_pairs
    fdupfree_mapped_pairs = 1.0 * dupfree_mapped_pairs / read_pairs
except OSError:
    exit("ERROR! Unable to read: {}".format(infile_MarkDuplicates))


try:
    # Picard InsertSizeMetrics
    if os.path.isfile(input_InsertSizeMetrics):
        with open(input_InsertSizeMetrics) as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].startswith("MEDIAN_INSERT_SIZE"):
                fragment_size = lines[i + 1].split()[0]
    # MACS2 xls
    elif os.path.isfile(infile_MACS2_xls):
        with open(infile_MACS2_xls) as f:
            lines = f.readlines()
        columns = list(filter(lambda x: x.startswith("# d"), lines))[0].strip()
        fragment_size = int(columns.split(" = ")[1])
    else:
        fragment_size = "NA"
except OSError:
    fragment_size = "NA"


# get peak qc from MACS2_peak_qc output #######################################
try:
    with open(infile_MACS2_qc_txt) as f:
        columns = list(map(lambda x: float(x), f.readlines()[1].split()))
        peak_count = int(columns[0])
        frip = round(float(columns[1]), 3)
        peak_genome_coverage = round(columns[2], 3)
except OSError:
    peak_count = "NA"
    frip = "NA"
    peak_genome_coverage = "NA"


# Output to stdout ###########################################################
print("{sample_name}\t{read_pairs}\t{mapped_pairs}\t{fmapped_pairs:.3f}\t{dup_mapped_pairs}\t{fdup_mapped_pairs:.3f}\t{dupfree_mapped_pairs}\t{fdupfree_mapped_pairs:.3f}\t{fhq_mapped_reads:.3f}\t{fmapped_singletons:.3f}\t{fragment_size}\t{peak_count}\t{frip}\t{peak_genome_coverage}".format(
    sample_name=os.path.basename(infile_AlignmentSummaryMetrics).split(".")[0],
    read_pairs=read_pairs,
    mapped_pairs=mapped_pairs,
    fmapped_pairs=fmapped_pairs,

    dup_mapped_pairs=dup_mapped_pairs,
    fdup_mapped_pairs=fdup_mapped_pairs,
    dupfree_mapped_pairs=dupfree_mapped_pairs,
    fdupfree_mapped_pairs=fdupfree_mapped_pairs,
    fhq_mapped_reads=fhq_mapped_reads,
    fmapped_singletons=fmapped_singletons,
    fragment_size=fragment_size,

    peak_count=peak_count,
    frip=frip,
    peak_genome_coverage=peak_genome_coverage,
))
