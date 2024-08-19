#!/usr/bin/env python
import argparse
import math

parser = argparse.ArgumentParser("Summarize a table of reads/feature/cell/UMI into tables of reads/feature/cell, UMIs/feature/cell and Poisson corrected UMIs/feature/cell")
parser.add_argument("--umiLength", type=int, default=4, help="UMI length")
parser.add_argument("input", help="Input count file (a tsv matrix of integer counts).")
parser.add_argument("ReadCounts", help="The output file name that will contain the number of reads per feature per cell")
parser.add_argument("UMICounts", help="The output file name that will contain the uncorrected number of UMIs per feature per cell")
parser.add_argument("CorrectedCounts", help="The output file name that will contain the corrected number of UMIs per feature per cell")
args = parser.parse_args()

maxUMIs = 4**args.umiLength
readCounts = dict()
UMICounts = dict()
lastGeneID = None

f = open(args.input)
header = None
for line in f:
    cols = line.split()
    if cols[0] == "GENEID":
        header = cols
        continue
    # Columns are now: Gene__chromosome UMI counts...
    if 'N' in cols[1]:
        continue
    if lastGeneID != cols[0]:
        seen = set()
        lastGeneID = cols[0]

    # Initialize to 0
    if cols[0] not in readCounts:
        readCounts[cols[0]] = [0] * (len(cols) - 2)
        UMICounts[cols[0]] = [0] * (len(cols) - 2)

    # Add the read counts
    readCounts[cols[0]] = [x + int(y) for x, y in zip(readCounts[cols[0]], cols[2:], strict=False)]

    for idx, cnt in enumerate(cols[2:]):
        if cnt == '0':
            continue
        if (cols[0], cols[1], idx) not in seen:
            UMICounts[cols[0]][idx] += 1
            seen.update((cols[0], cols[1], idx))
f.close()

ReadCountsFile = open(args.ReadCounts, "w")
UMICountsFile = open(args.UMICounts, "w")
CorrectedCountsFile = open(args.CorrectedCounts, "w")

# Write the output
for f in [ReadCountsFile, UMICountsFile, CorrectedCountsFile]:
    f.write("GENEID\t{}\n".format("\t".join(header[2:])))

sortedFeatures = sorted(readCounts.keys())
for feature in sortedFeatures:
    ReadCountsFile.write("{}\t{}\n".format(feature, "\t".join(["{}".format(x) for x in readCounts[feature]])))
    UMICountsFile.write("{}\t{}\n".format(feature, "\t".join(["{}".format(x) for x in UMICounts[feature]])))
    CorrectedCountsFile.write("{}".format(feature))
    for cnt in UMICounts[feature]:
        if (cnt >= maxUMIs):
            maxUMIs -= 0.5
        CorrectedCountsFile.write("\t{}".format(abs(-math.log(1.0 - float(cnt) / float(maxUMIs)) * maxUMIs)))  # The abs() prevents -0.0
    CorrectedCountsFile.write("\n")
for f in [ReadCountsFile, UMICountsFile, CorrectedCountsFile]:
    f.close()
