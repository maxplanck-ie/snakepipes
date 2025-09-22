#!/usr/bin/env python
'''Convert a bedGraph file from MethylDackel to two bigWig files'''
import pyBigWig
import sys
import os

fai = open(sys.argv[2])
hdr = []
for line in fai:
    cols = line.strip().split()
    hdr.append((cols[0], int(cols[1])))
fai.close()

bg = open(sys.argv[1])
bwM = pyBigWig.open(sys.argv[3], "w")
bwC = pyBigWig.open(sys.argv[4], "w")

bwM.addHeader(hdr)
bwC.addHeader(hdr)

chrom = []
starts = []
ends = []
valuesM = []
valuesC = []
for line in bg:
    if line.startswith("track "):
        continue
    cols = line.strip().split()
    if len(chrom) and chrom[0] != cols[0]:
        bwM.addEntries(chrom, starts, ends=ends, values=valuesM)
        bwC.addEntries(chrom, starts, ends=ends, values=valuesC)
        chrom = []
        starts = []
        ends = []
        valuesM = []
        valuesC = []
    chrom.append(cols[0])
    starts.append(int(cols[1]))
    ends.append(int(cols[2]))
    valuesM.append(float(cols[3]))
    valuesC.append(float(cols[4]) + float(cols[5]))

if chrom:
    bwM.addEntries(chrom, starts, ends=ends, values=valuesM)
    bwC.addEntries(chrom, starts, ends=ends, values=valuesC)

bg.close()
bwM.close()
bwC.close()
