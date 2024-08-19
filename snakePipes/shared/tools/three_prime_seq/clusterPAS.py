#!/usr/bin/env python
import argparse



def cluster_pas(args):
    """
    Clusters the PAS and defines position with highest signal.
    """

    window = int(args.windowSize)
    minReads = int(args.minReads)

    clusters = {}

    with open(args.infile, 'r') as In:
        #line = In.readline()
        line = In.readline()
        while line:
            entry = line.strip("\n").split('\t')
            chrom = entry[0]
            pas = int(entry[1])
            strand = entry[2]
            count = float(entry[3])
            gene = entry[4]
            annotation = entry[7]

            if count > minReads:
                if gene not in clusters:
                    clusters[gene] = {}
                    #cluster start, cluster end, max count, gene, annotation, summit
                    clusters[gene][pas] = [pas, pas + 1, count, gene, annotation, pas, strand, chrom]
                else:
                    newCluster = True
                    for c in clusters[gene]:
                        if pas > clusters[gene][c][0] - window and pas < clusters[gene][c][1] + window:
                            clusters[gene][c][0] = min(pas, clusters[gene][c][0])
                            clusters[gene][c][1] = max(pas + 1, clusters[gene][c][1])
                            clusters[gene][c][2] += count
                            if count > clusters[gene][c][2]:
                                #clusters[gene][c][2] = count
                                clusters[gene][c][4] = annotation
                                clusters[gene][c][5] = pas
                            newCluster = False
                            break

                    #not close to any existing cluster => new cluster
                    if newCluster:
                        clusters[gene][pas] = [pas, pas + 1, count, gene, annotation, pas, strand, chrom]

            line = In.readline()

    return clusters


def parseArgs():
    parser = argparse.ArgumentParser(description="Clusters poly(A) sites that are supported by at least m reads within a given distance w. M can be changed with --minReads, w with --windowSize.")
    parser.add_argument("infile", help="Input file with annotated PAS.")
    parser.add_argument("output", help="Output file with clusters of PAS.")
    parser.add_argument("--windowSize", help="Maximal distance between PAS to be clustered (default: 15).", default=15, type=int)
    parser.add_argument("--minReads", help="Number of reads that define a PAS (default: 5).", default=5, type=int)
    return parser


# Input: chrom\tpos\tstrand\tcount\tgene,list\ttranscript,list\ttranscript_positions,list\ttag (UTR, exon, extension)
# Output: Chromosome\tStart\tEnd\tGene\tCounts\tStrand\tAnnotation\tSummit\n
def main():
    args = parseArgs().parse_args()

    clusters = cluster_pas(args)

    o = open(args.output, "w")
    o.write("Chromosome\tStart\tEnd\tGene\tCounts\tStrand\tAnnotation\tSummit\n")
    for chrom in clusters:
        for clust in clusters[chrom]:
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(clusters[chrom][clust][7], clusters[chrom][clust][0], clusters[chrom][clust][1], clusters[chrom][clust][3], clusters[chrom][clust][2], clusters[chrom][clust][6],clusters[chrom][clust][4],clusters[chrom][clust][5]))
    o.close()


if __name__ == "__main__":
    main()
