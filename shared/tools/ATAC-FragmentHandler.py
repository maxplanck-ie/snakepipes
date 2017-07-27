import sys

class AtacFragment:
    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

        self.size = end-start

    def shiftCutsite(self):
        self.start -= 4
        self.end += 5

    def formatBed(self):
        return('\t'.join([self.chrom, str(self.start), str(self.end), '', '', self.strand]))


class ChromosomeBoundaries:
    def __init__(self, chromosomes, sizes):
        self.chromMap = {}
        for i in range(len(chromosomes)):
            self.chromMap[chromosomes[i]] = sizes[i]

    def isInBoundary(self, fragment):
        offset_min = 0
        offset_max = 0

        maxend = self.chromMap[fragment.chrom]
        if fragment.start < 0:
            offset_min = 1
        if fragment.end > maxend:
            offset_max = 1
        return(offset_min, offset_max)

def main():
    fastaindex = sys.argv[2]
    with open(fastaindex,'r') as fai:
        tmp = fai.readlines()
        tmp = [x.split('\t') for x in tmp]
        chroms = [x[0] for x in tmp]
        sizes = [int(x[2]) for x in tmp]
    chromHandler = ChromosomeBoundaries(chroms, sizes)

    for line in sys.stdin:
        tmp = line.split('\t')
        chrom, start, end, strand = [tmp[i] for i in (0,1,5,8)]
        frag = AtacFragment(chrom, int(start), int(end), strand)
        lb, ub = chromHandler.isInBoundary(frag)
        if (lb):
            print('lb')
        if (ub):
            print('ub')
        print(frag)

if __name__=='__main__':
    try:
        main()
    except IOError:
        sys.stderr.write("IOError" + '\n')
    except:
        sys.stderr.write("Unknown error" + '\n')
        raise
