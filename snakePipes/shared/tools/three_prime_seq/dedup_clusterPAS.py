#!/usr/bin/env python

import pandas as pd
import sys

"""
Read in post-processed output of clusterPAS from STDIN, add header,
deduplicate the Gene column (4th), write to STDOUT.
"""


HEADERS = [
    "Chromosome",
    "Start",
    "End",
    "Gene",
    "Counts",
    "Strand",
    "Annotation",
    "Summit",
]


def dedup(sub_df):
    new_ids = ["_".join([str(gene), str(i)]) for
               (i, gene) in enumerate(sub_df['Gene'])]
    sub_df['Gene'] = new_ids
    return sub_df


def main():
    df = pd.read_table(sys.stdin, index_col=None, header=0)
    df = (
        df.groupby("Gene", group_keys=False)
        .apply(dedup)
        .reset_index(drop=True)
    )
    df.to_csv(sys.stdout, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
