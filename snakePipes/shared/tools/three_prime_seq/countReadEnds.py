#!/usr/bin/env python
import argparse
import pyBigWig
import numpy as np


def count_readEnds(args):
    counts = {}
    bwFw = pyBigWig.open(args.forward)
    bwRev = pyBigWig.open(args.reverse)
    for line in open(args.bedfile):
        cols = line.strip().split()
        #print(line)
        #strand specific signal
        if cols[5] == "+":
            vals = bwFw.values(cols[0], int(cols[1]), int(cols[2]))
        else:
            vals = bwRev.values(cols[0], int(cols[1]), int(cols[2]))
        countSum = np.nansum(vals)
        #save in counts
        if cols[3] in counts:
            c = counts[cols[3]][4]
            counts[cols[3]] = [cols[0], cols[1], cols[2], cols[3], c + countSum, cols[5]]
        else:
            counts[cols[3]] = [cols[0], cols[1], cols[2], cols[3], countSum, cols[5]]
    bwFw.close()
    bwRev.close()
    return counts



def parseArgs():
    parser = argparse.ArgumentParser(description="Counts number of 5' ends in each region of a BED-file.")
    parser.add_argument("forward", help="Strand-specific bigwig file (forward strand) with 5' signal of _R1.")
    parser.add_argument("reverse", help="Strand-specific bigwig file (reverse strand) with 5' signal of _R1.")
    parser.add_argument("bedfile", help="Input bed file.")
    parser.add_argument("output", help="Output file number of 5' read counts for each region.")
    return parser



def main():
    args = parseArgs().parse_args()

    readCounts = count_readEnds(args)

    o = open(args.output, "w")
    o.write("Chromosome\tStart\tEnd\tGene\tCounts\tStrand\n")
    for region in readCounts:
        o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(readCounts[region][0], readCounts[region][1], readCounts[region][2], readCounts[region][3], int(readCounts[region][4]), readCounts[region][5]))
    o.close()


if __name__ == "__main__":
    main()




#paste Sample1.binnedTandemUTRs.counts.txt Sample2.binnedTandemUTRs.counts.txt Sample3.binnedTandemUTRs.counts.txt Sample4.binnedTandemUTRs.counts.txt Sample5.binnedTandemUTRs.counts.txt Sample6.binnedTandemUTRs.counts.txt Sample7.binnedTandemUTRs.counts.txt Sample8.binnedTandemUTRs.counts.txt Sample9.binnedTandemUTRs.counts.txt Sample10.binnedTandemUTRs.counts.txt Sample11.binnedTandemUTRs.counts.txt Sample12.binnedTandemUTRs.counts.txt | awk '{OFS="\t"; print $1,$2,$3,$4,".",$6,$5,$11,$17,$23,$29,$35,$41,$47,$53,$59,$65,$71}' > NewPipeline_CountReadEnds_EndToEnd_allSamples.txt
