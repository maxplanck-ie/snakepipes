#!/usr/bin/env python3


# single-end ONLY!


# Usage: python sample_qc_report.py sample_name.alignment_summary_metrics.txt sample_name.mark_duplicates_metrics.txt sample_name.filtered.BAM_peaks.xls sample_name.filtered.BAM_peaks.qc.txt

import os
import sys


# Input files ################################################################
infile_AlignmentSummaryMetrics = sys.argv[1]    # sample_name.alignment_summary_metrics.txt
infile_MarkDuplicates = sys.argv[2]             # sample_name.mark_duplicates_metrics.txt
try:
    infile_MACS2_xls = sys.argv[3]              # sample_name.filtered.BAM_peaks.xls (optional)
    infile_MACS2_qc_txt = sys.argv[4]           # ample_name.filtered.BAM_peaks.qc.txt (optional)
except IndexError:
    pass


# read Picard AlignmentSummaryMetrics results ################################
try:
    with open(infile_AlignmentSummaryMetrics) as f:
        lines = f.readlines()
    columns = list(filter(lambda x: x.startswith("UNPAIRED"), lines))[0].strip().split("\t")

    # PF_READS: The number of PF reads where PF is defined as passing Illumina's filter.
    total_reads = int(columns[2])

    # PF_READS_ALIGNED: The number of PF reads that were aligned to the reference sequence.
    mapped_reads = int(columns[5])

    # PCT_PF_READS_ALIGNED: The percentage of PF reads that aligned to the reference sequence. PF_READS_ALIGNED / PF_READS
    fmapped_reads = float(columns[6])

    # fraction of high-quality mapped reads (MAPQ >=20)
    fhq_mapped_reads = 1.0 * int(columns[8]) / total_reads
except OSError:
    exit("ERROR! Unable to read: {}\n".format(infile_AlignmentSummaryMetrics))


# read Picard MarkDuplicates results #########################################
try:
    with open(infile_MarkDuplicates) as f:
        lines = f.readlines()
    columns = list(filter(lambda x: x.startswith("Unknown Library"), lines))[0].strip().split("\t")

    # UNPAIRED_READ_DUPLICATES: The number of fragments that were marked as duplicates.
    dup_mapped_reads = int(columns[4])

    # fraction of UNPAIRED_READ_DUPLICATES
    fdup_mapped_reads = 1.0 * dup_mapped_reads / mapped_reads

    # duplication free mapped reads
    dupfree_mapped_reads = mapped_reads - dup_mapped_reads

    # fraction of duplication free mapped reads
    fdupfree_mapped_reads = 1.0 * dupfree_mapped_reads / total_reads
except OSError:
    exit("ERROR! Unable to read: {}".format(infile_MarkDuplicates))


# get fragment size from Picard CollectInsertSizeMetrics/MACS2 output ########
try:
    with open(infile_MACS2_xls) as f:
        lines = f.readlines()
    columns = list(filter(lambda x: x.startswith("# d"), lines))[0].strip()
    fragment_size = int(columns.split(" = ")[1])
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
print("{sample_name}\t{total_reads}\t{mapped_reads}\t{fmapped_reads:.3f}\t{dup_mapped_reads}\t{fdup_mapped_reads}\t{dupfree_mapped_reads}\t{fdupfree_mapped_reads}\t{fhq_mapped_reads:.3f}\t{fragment_size}\t{peak_count}\t{frip}\t{peak_genome_coverage}".format(
    sample_name=os.path.basename(infile_AlignmentSummaryMetrics).split(".alignment_summary_metrics.txt")[0],
    total_reads=total_reads,
    mapped_reads=mapped_reads,
    fmapped_reads=fmapped_reads,
    dup_mapped_reads=dup_mapped_reads,
    fdup_mapped_reads=fdup_mapped_reads,
    dupfree_mapped_reads=dupfree_mapped_reads,
    fdupfree_mapped_reads=fdupfree_mapped_reads,
    fhq_mapped_reads=fhq_mapped_reads,
    fragment_size=fragment_size,
    peak_count=peak_count,
    frip=frip,
    peak_genome_coverage=peak_genome_coverage))
