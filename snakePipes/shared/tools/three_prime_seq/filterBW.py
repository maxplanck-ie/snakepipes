#!/usr/bin/env python
from deeptoolsintervals import GTF
import pyBigWig
import argparse
import os

def filterOverlaps(chrom, interval, BED):
    """
    Remove portions of an interval overlapping a BED file, returning a list of
    (start, end, value) tuples.

    Also returns the number of bases filtered
    """
    s, e, v = interval
    olaps = BED.findOverlaps(chrom, s, e)
    if not olaps or len(olaps) == 0:
        return 0, [interval]
    blacklisted = set()
    for o in olaps:
        for pos in range(o[0], o[1]):
            blacklisted.add(pos)

    finalBases = [x for x in range(s, e) if x not in blacklisted]
    filteredBases = (e - s) - len(finalBases)

    if not len(finalBases):
        return filteredBases, []

    # Convert that to intervals so things are better compressed
    lastStart = finalBases[0]
    lastEnd = lastStart + 1
    finalIntervals = []
    for base in finalBases[1:]:
        if base == lastEnd:
            lastEnd = base
        else:
            finalIntervals.append((lastStart, lastEnd, v))
            lastStart, lastEnd = base, base + 1
    finalIntervals.append((lastStart, lastEnd, v))
    return filteredBases, finalIntervals


parser = argparse.ArgumentParser(description="Blacklist signal from a bigWig file")
parser.add_argument("bigWig", help="Input bigWig file")
parser.add_argument("BED", help="Input BED file")
parser.add_argument("output", help="Output bigWig file")
args = parser.parse_args()

bw = pyBigWig.open(args.bigWig)
BED = GTF(args.BED)
o = pyBigWig.open(args.output, "w")
filteredBases = 0
totalBases = bw.header()['nBasesCovered']

# write the header
hdr = [(k, v) for k, v in bw.chroms().items()]
o.addHeader(hdr)

for chrom, _ in hdr:
    ints = bw.intervals(chrom)
    filteredStarts = []
    filteredEnds = []
    filteredValues = []
    if not ints:
        continue
    for i in ints:
        fbases, flist = filterOverlaps(chrom, i, BED)
        filteredBases += fbases
        for s, e, v in flist:
            filteredStarts.append(s)
            filteredEnds.append(e)
            filteredValues.append(v)
    if len(filteredStarts):
        o.addEntries([chrom] * len(filteredStarts), filteredStarts, ends=filteredEnds, values=filteredValues)
bw.close()
o.close()

print("Sample\tFilteredBases\tTotalBases")
print("{}\t{}\t{}".format(os.path.basename(args.bigWig)[:-3], filteredBases, totalBases))
