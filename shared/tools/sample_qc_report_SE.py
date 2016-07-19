#!/usr/bin/env python


## single-end ONLY!


## Usage: python sample_qc_report.py sample_name.alignment_summary_metrics.txt sample_name.mark_duplicates_metrics.txt sample_name.filtered.BAM_peaks.xls sample_name.filtered.BAM_peaks.qc.txt

import os
import sys


## Input files ################################################################
infile_AlignmentSummaryMetrics = sys.argv[1]    # sample_name.alignment_summary_metrics.txt
infile_MarkDuplicates = sys.argv[2]             # sample_name.mark_duplicates_metrics.txt
try:
    infile_MACS2_xls = sys.argv[3]              # sample_name.filtered.BAM_peaks.xls (optional)
    infile_MACS2_qc_txt = sys.argv[4]           # ample_name.filtered.BAM_peaks.qc.txt (optional)
except:
    pass


## read Picard AlignmentSummaryMetrics results ################################
try:
    with open(infile_AlignmentSummaryMetrics) as f:
        lines = f.readlines()
    columns = filter(lambda x: x.startswith("UNPAIRED"), lines)[0].strip().split("\t")
    ##print columns

    ## PF_READS: The number of PF reads where PF is defined as passing Illumina's filter.
    total_reads = int(columns[2])
    ##print "total_reads:", total_reads

    ## PF_READS_ALIGNED: The number of PF reads that were aligned to the reference sequence.
    mapped_reads = int(columns[5])
    ##print "mapped_reads:", mapped_reads

    ## PCT_PF_READS_ALIGNED: The percentage of PF reads that aligned to the reference sequence. PF_READS_ALIGNED / PF_READS
    fmapped_reads = float(columns[6])
    ##print "fmapped_reads:", fmapped_reads

    ## fraction of high-quality mapped reads (MAPQ >=20)
    fhq_mapped_reads = 1.0 * int(columns[8]) / total_reads
    ##print "fhq_mapped_reads:", fhq_mapped_reads
except:
    print "ERROR! Unable to read:", infile_AlignmentSummaryMetrics
    exit(1)


## read Picard MarkDuplicates results #########################################
try:
    with open(infile_MarkDuplicates) as f:
        lines = f.readlines()
    columns = filter(lambda x: x.startswith("Unknown Library"), lines)[0].strip().split("\t")
    ##print columns

    ## UNPAIRED_READ_DUPLICATES: The number of fragments that were marked as duplicates.
    dup_mapped_reads = int(columns[4])
    ##print "dup_mapped_reads:", dup_mapped_reads

    ## fraction of UNPAIRED_READ_DUPLICATES
    fdup_mapped_reads = 1.0 * dup_mapped_reads / mapped_reads
    ##print "fdup_mapped_reads:", fdup_mapped_reads

    ## duplication free mapped reads
    dupfree_mapped_reads = mapped_reads - dup_mapped_reads
    ##print "dupfree_mapped_reads:", dupfree_mapped_reads

    ## fraction of duplication free mapped reads
    fdupfree_mapped_reads = 1.0 * dupfree_mapped_reads / total_reads
    ##print "fdupfree_mapped_reads:", fdupfree_mapped_reads
except:
    print "ERROR! Unable to read:", infile_MarkDuplicates
    exit(1)


## get fragment size from Picard CollectInsertSizeMetrics/MACS2 output ########
try:
    with open(infile_MACS2_xls) as f:
        lines = f.readlines()
    columns = filter(lambda x: x.startswith("# d"), lines)[0].strip()
    fragment_size = int(columns.split(" = ")[1])
except:
    fragment_size = "NA"
##print "fragment_size:", fragment_size


# get peak qc from MACS2_peak_qc output #######################################
try:
    with open(infile_MACS2_qc_txt) as f:
        peak_count, frip, peak_genome_coverage = f.readlines()[1].split()
except:
    peak_count = "NA"
    frip = "NA"
    peak_genome_coverage = "NA"
##print "peak_count:", peak_count
##print "frip:", frip
##print "peak_genome_coverage:", peak_genome_coverage


## Output to stdout ###########################################################
print "{sample_name}\t{total_reads}\t{mapped_reads}\t{fmapped_reads:.3f}\t{dup_mapped_reads}\t{fdup_mapped_reads}\t{dupfree_mapped_reads}\t{fdupfree_mapped_reads}\t{fhq_mapped_reads:.3f}\t{fragment_size}\t{peak_count}\t{frip}\t{peak_genome_coverage}".format(
    sample_name=os.path.basename(infile_AlignmentSummaryMetrics).split(".")[0],
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
    peak_genome_coverage=peak_genome_coverage)
