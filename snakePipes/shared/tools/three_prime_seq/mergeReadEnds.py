#!/usr/bin/env python

import pandas as pd
import sys
import argparse

"""
Munge multiple outputs of countReadEnds into a single counts.tsv
"""

def munge(infiles, samples):
    final_df = list()
    for fn, sample in zip(infiles, samples, strict=False):
        df = pd.read_table(fn, header=0, index_col=None)
        df['Sample'] = sample
        df = df[['Gene', 'Counts', 'Sample']]
        final_df.append(df)
    final_df = pd.concat(final_df)
    final_df = pd.pivot_table(final_df, index="Gene", values="Counts", columns="Sample")
    final_df.index.name = ''
    return final_df



def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--infiles", type=argparse.FileType("r"), nargs="+")
    parser.add_argument('-s', "--samples", type=str, nargs='+')
    parser.add_argument("-o", "--outfile", type=argparse.FileType("w"), default=sys.stdout)
    args = parser.parse_args(argv)
    if len(args.infiles) != len(args.samples):
        raise Exception('Error: infiles and samples length mismatch.')
    df = munge(args.infiles, args.samples)
    df.to_csv(args.outfile, sep="\t")


if __name__ == "__main__":
    main(sys.argv[1:])
