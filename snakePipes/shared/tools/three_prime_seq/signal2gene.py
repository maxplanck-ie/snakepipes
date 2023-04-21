#!/usr/bin/env python
import argparse
import os.path
from deeptoolsintervals import GTF
from deeptoolsintervals.parse import tree, openPossiblyCompressed, findRandomLabel
from deeptools.mapReduce import mapReduce
import pyBigWig
import csv

class extendedGTF(GTF):
    def __init__(self, fname, extend=500):
        self.fname = []
        self.filename = ""
        self.chroms = []
        self.exons = []
        self.labels = []
        self.transcriptIDduplicated = []
        self.tree = tree.initTree()
        self.labelIdx = 0
        self.transcript_id_designator = "transcript_id"
        self.gene_id_designator = "gene_id"
        self.exonID = "exon"
        self.transcriptID = "transcript"
        self.keepExons = True
        self.defaultGroup = None
        self.verbose = False
        self.extend = extend
        self.transcripID2Gene = dict()
        self.labelIdx = 0

        self.filename = fname
        fp = openPossiblyCompressed(fname)
        line, labelColumn = self.firstNonComment(fp)
        assert(line)
        line = line.strip()

        self.parseGTF(fp, line)
        # vine -> tree
        self.tree.finish()

    def parseGID(self, attribs):
        """
        Return the gene_id and transcript_id from the attributes
        """
        s = next(csv.reader([attribs], delimiter=' '))
        tid = s[s.index(self.transcript_id_designator) + 1].rstrip(";")
        gid = s[s.index(self.gene_id_designator) + 1].rstrip(";")
        return gid, tid

    def parseGTF(self, fp, line):
        """
        The only difference is that this adds args.extension onto the appropriate end
        and creates a transcriptID -> geneID map
        """
        file_label = findRandomLabel(self.labels, os.path.basename(self.filename))

        # Handle the first line
        cols = line.split("\t")
        if cols[2].lower() == self.transcriptID.lower():
            self.parseGTFtranscript(cols, file_label)
        elif cols[2].lower() == self.exonID.lower():
            self.parseGTFexon(cols)

        # Handle the remaining lines
        for line in fp:
            if not isinstance(line, str):
                line = line.decode('ascii')
            if not line.startswith('#'):
                cols = line.split("\t")
                if len(cols) == 0:
                    continue

                if cols[2].lower() == self.transcriptID.lower():
                    if cols[6] == "-":
                        cols[3] = "{}".format(max(1, int(cols[3]) - self.extend))
                    else:
                        cols[4] = "{}".format(int(cols[4]) + self.extend)
                    self.parseGTFtranscript(cols, file_label)

                    # Add transcript_id -> gene_id mapping
                    gid, tid = self.parseGID(cols[8])
                    self.transcripID2Gene[tid] = gid
                elif cols[2].lower() == self.exonID.lower() and self.keepExons is True:
                    self.parseGTFexon(cols)

        # Reset self.labelIdx
        self.labelIdx = len(self.labels)

def o2attribs(chrom, pos, o):
    tid = o[2]
    geneID = exonsTree.transcripID2Gene[tid]
    tag = None

    # Does this actually overlap an exon?
    for exon in o[4]:
        if pos >= exon[0] and pos < exon[1]:
            tag = "exonic"
            break
    if not tag:
        if o[5] == "-":
            if pos < o[4][0][0]:
                tag = "downstream"
                transPos = "+{}".format(o[4][0][0] - pos)
        else:
            if pos >= o[4][-1][1]:
                tag = "downstream"
                transPos = "+{}".format(pos - o[4][-1][1])

    # Add UTR and CDS tags
    if tag == "exonic":
        # CDS, 3'UTR, or 5'UTR
        cdss = CDS.exons[0][tid]
        if len(cdss) > 0:
            for cds in cdss:
                if pos >= cds[0] and pos < cds[1]:
                    tag = "CDS"
                    break
            if tag != "CDS":
                if o[5] == "-":
                    if pos < cdss[0][0]:
                        tag = "3'UTR"
                    else:
                        tag = "5'UTR"
                else:
                    if pos >= cdss[-1][1]:
                        tag = "3'UTR"
                    else:
                        tag = "5'UTR"

    if not tag:
        return None, None, None, None

    # Calculate position in transcript coordinates
    if tag != "downstream":
        transPos = 0
        if o[5] == "-":
            for exon in o[4][::-1]:
                if pos > exon[1]:
                    break
                if pos < exon[0]:
                    transPos += exon[1] - exon[0]
                else:
                    transPos += exon[1] - pos
                    break
        else:
            for exon in o[4]:
                if pos < exon[0]:
                    break
                if pos > exon[1]:
                    transPos += exon[1] - exon[0]
                else:
                    transPos += pos - exon[0]
                    break

    return geneID, tid, transPos, tag


def mapStrand(bw, chrom, start, end, strand, args):
    res = []
    intervals = bw.intervals(chrom, start, end)
    if not intervals:
        return res
    for interval in intervals:
        if interval[0] < start:
            continue
        for pos in range(interval[0], interval[1]):
            o = exonsTree.findOverlaps(chrom, pos, pos + 1, strand=strand, strandType=3, includeStrand=True)
            genes = set()
            transcripts = []
            transPos = []
            tags = []
            for overlap in o:
                geneID, transcriptID, tPos, tag = o2attribs(chrom, pos, overlap)
                if not geneID:
                    continue  # intronic
                genes.add(geneID)
                transcripts.append(transcriptID)
                transPos.append("{}".format(tPos))
                tags.append(tag)

            if len(genes):
                res.append([chrom, pos, strand, interval[2], genes, transcripts, transPos, tags])

    return res


def mapValuesWrapper(args):
    return mapValues_worker(*args)


def mapValues_worker(chrom, start, end, args):
    """
    Actually extract signal and annotate them.
    """
    bwF = pyBigWig.open(args.forward_bigWig)
    bwR = pyBigWig.open(args.reverse_bigWig)

    res = mapStrand(bwF, chrom, start, end, "+", args)
    res.extend(mapStrand(bwR, chrom, start, end, "-", args))
    bwF.close()
    bwR.close()
    return res


def parseArgs():
    parser = argparse.ArgumentParser(description="Associate signal with position on a gene or genes.")
    parser.add_argument("forward_bigWig", help="Input BigWig file for the forward strand.")
    parser.add_argument("reverse_bigWig", help="Input BigWig file for the reverse strand.")
    parser.add_argument("GTF", help="Input GTF file.")
    parser.add_argument("output", help="Output file.")
    parser.add_argument("--extend", help="Number of bases to extend each gene in the 3' direction (default: 100).", default=100, type=int)
    parser.add_argument("--threads", help="Number of threads (default: 1).", default=1, type=int)
    return parser

# Load GTF as normal
# Use transcript -> gene mapping
# Output: chrom\tpos\tstrand\tcount\tgene,list\ttranscript,list\ttranscript_positions,list\ttag (UTR, exon, extension)
def main():
    args = parseArgs().parse_args()

    bwF = pyBigWig.open(args.forward_bigWig)
    chromList = [(k, v) for k, v in bwF.chroms().items()]
    bwF.close()

    global exonsTree, CDS
    CDS = GTF(args.GTF, keepExons=True, exonID="CDS")
    exonsTree = extendedGTF(args.GTF, extend=args.extend)  # Adds exonsTree.transcriptID2Gene and 3'UTR extension
    res = mapReduce([args], mapValuesWrapper, chromList, genomeChunkLength=1000000, numberOfProcessors=args.threads)

    o = open(args.output, "w")
    o.write("Chromosome\tPosition\tStrand\tCounts\tGenes\tTranscripts\tTranscriptPositions\tTags\n")
    for block in res:
        for r in block:
            line = [r[0], r[1], r[2], r[3], ",".join(list(r[4])), ",".join(r[5]), ",".join(r[6]), ",".join(r[7])]
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*line))
    o.close()

if __name__ == "__main__":
    main()
