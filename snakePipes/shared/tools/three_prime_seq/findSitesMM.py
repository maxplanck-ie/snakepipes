#!/usr/bin/env python
# This can take up to ~20 minutes and use up to ~2GB RAM for mammals
import argparse
import py2bit
from deeptoolsintervals import GTF, tree
from deeptoolsintervals.parse import openPossiblyCompressed, parseExonBounds, findRandomLabel
import sys

parser = argparse.ArgumentParser(description="Generate a blacklist file of polyX stretches of a given minimum length not within a specified distance of a TES")
parser.add_argument("--output", "-o", help="Output file", required=True)
parser.add_argument("--tb", help="2bit file", required=True)
parser.add_argument("--bed", help="BED file containing transcripts", required=True)
parser.add_argument("--windowLength", help="Minimum length of downstream window (default: 10)", type=int, default=10)
parser.add_argument("--minLength", help="Minimum length of a polyX stretch (default: 6)", type=int, default=6)
parser.add_argument("--base", help="The base to check for (A or T)", required=True, choices=['A', 'T'])
parser.add_argument("--percBase", help="Pecentage of base (A or T) required to be present in downstream window (default: 0.7)", type=float, default=0.7)
parser.add_argument("--minDistance", help="Minimum distance from a TES to not exclude a site (default: 100)", type=int, default=100)
parser.add_argument("--extend", help="Number of bases to extend regions (default: 5)", type=int, default=5)
args = parser.parse_args()

tb = py2bit.open(args.tb)
o = open(args.output, "w")

class TES(GTF):
    def parseBEDcore(self, line, ncols):
        strand = 3
        cols = line.split("\t")
        name = "{0}:{1}-{2}".format(cols[0], cols[1], cols[2])

        if int(cols[1]) < 0:
            cols[1] = 0

        if int(cols[1]) >= int(cols[2]):
            sys.stderr.write("Warning: {0}:{1}-{2} is an invalid BED interval! Ignoring it.\n".format(cols[0], cols[1], cols[2]))
            return

        # BED6/BED12: set name and strand
        score = '.'
        if ncols > 3:
            name = cols[3]
            if cols[5] == '+':
                strand = 0
            elif cols[5] == '-':
                strand = 1
            score = cols[4]

        # filter by strand
        if strand != 3:
            if self.strand == "+" and strand == 1:
                return
            elif self.strand == "-" and strand == 0:
                return

        # Ensure that the name is unique
        name = findRandomLabel(self.exons[self.labelIdx], name)

        assert(len(cols) == 12)
        exons = parseExonBounds(int(cols[1]), int(cols[2]), int(cols[9]), cols[10], cols[11])

        # Extend by strand around the TES
        exonsFinal = []
        lenLeft = 2 * self.minDistance + 1
        if strand == 3 or (self.strand == "+" and strand == 0):
            exons[-1] = (exons[-1][0], exons[-1][1] + self.minDistance)
            for exon in exons[::-1]:
                exonLen = exon[1] - exon[0]
                if exonLen <= lenLeft:
                    lenLeft -= exonLen
                    exonsFinal.insert(0, exon)
                else:
                    exonsFinal.insert(0, (exon[1] - lenLeft, exon[1]))
                    break
                if lenLeft <= 0:
                    break
        elif self.strand == "-" and strand == 1:
            _ = exons[0][0] -self.minDistance
            _ = max(0, _)
            exons[0] = (_, exons[0][1])
            for exon in exons:
                exonLen = exon[1] - exon[0]
                if exonLen <= lenLeft:
                    lenLeft -= exonLen
                    exonsFinal.append(exon)
                else:
                    exonsFinal.append((exon[0], exon[0] + lenLeft))
                    break
                if lenLeft <= 0:
                    break
        if len(exonsFinal) == 0:
            return

        self.tree.addEntry(self.mungeChromosome(cols[0]), exonsFinal[0][0], exonsFinal[-1][1], name, strand, self.labelIdx, score)
        self.exons[self.labelIdx][name] = exonsFinal


    def __init__(self, fname, minDistance=100, strand="+"):
        self.fname = [fname]
        self.filename = fname
        self.filename = ""
        self.chroms = []
        self.exons = []
        self.labels = []
        self.transcriptIDduplicated = []
        self.tree = tree.initTree()
        self.labelIdx = 0
        self.keepExons = True
        self.defaultGroup = None
        self.verbose = False
        self.minDistance = minDistance
        self.strand = strand

        fp = openPossiblyCompressed(fname)
        line, labelColumn = self.firstNonComment(fp)
        assert(line)  # This will only fail on empty files
        line = line.strip()

        self.ftype = self.inferType(fp, line, labelColumn)
        self.parseBED(fp, line, 12, labelColumn)
        fp.close()

        # Sanity check
        if self.tree.countEntries() == 0:
            raise RuntimeError("None of the input BED/GTF files had valid regions")

        # vine -> tree
        self.tree.finish()


def processLast(last, chrom, idx, idx2, o, bed):
    if args.base == "A":
        s = max(0, idx - args.extend)
        e = idx2
    else:
        s = idx
        e = idx2 + args.extend
    for overlaps in bed.findOverlaps(chrom, s, e):
        for exon in overlaps[4]:
            if exon[0] < s and exon[1] > s:
                return
            if exon[0] >= s and exon[0] < e:
                return
    if not last[0]:
        last[0] = chrom
        last[1] = s
        last[2] = e
    else:
        if last[0] == chrom and s <= last[2]:
            last[2] = e
        else:
            o.write("{}\t{}\t{}\n".format(*last))
            last[0] = chrom
            last[1] = s
            last[2] = e




if args.base == "A":
    bed = TES(args.bed, minDistance=args.minDistance)
else:
    bed = TES(args.bed, minDistance=args.minDistance, strand="-")

last =  [None, None, None]

for chrom, chromLength in tb.chroms().items():
    s = tb.sequence(chrom)

    idx = 0
    idx2 = 0
    while idx < chromLength - args.windowLength:
       idx2 = idx + args.windowLength
       if s[idx:idx2].count(args.base)/args.windowLength >= args.percBase  or s[idx:idx+args.minLength+1].count(''.join(args.base*args.minLength)) > 0:
           if args.base == "A":
               processLast(last, chrom, idx, idx + 1, o, bed)
           else:
               processLast(last, chrom, idx2 -1, idx2, o, bed)
       idx += 1



if last[0] is not None:
    o.write("{}\t{}\t{}\n".format(*last))
o.close()
tb.close()
