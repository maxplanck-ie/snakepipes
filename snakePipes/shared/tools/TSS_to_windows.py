import pysam
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-inf", "--input", dest="inf", help="input bed file", metavar="FILE")
parser.add_argument("-outf", "--output", dest="outf", help="output bed file", metavar="FILE")
parser.add_argument("-bam", "--bam", dest="bam", help="bam file", metavar="FILE")
parser.add_argument("-size", "--size", dest="size", help="flank width", type=int)


args = parser.parse_args()

inf = args.inf
outf = args.outf
bamf = args.bam
size = args.size

bam = pysam.AlignmentFile(bamf)
chroms_sizes = dict(zip(bam.references, bam.lengths, strict=False))

with open(inf) as f, open(outf, 'w') as of:
    for line in f:
        linesplit = line.split('\t')
        chr = linesplit[0]
        if chr in chroms_sizes.keys():
            start = max(1, int(linesplit[1]) - size)
            end = min(int(chroms_sizes[chr]), int(linesplit[1]) + size)
            of.write("{}\t{}\t{}\n".format(chr, str(start), str(end)))
