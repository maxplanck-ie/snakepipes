#!/usr/bin/env python3

import sys
import gzip

def merge_pairs(cis1_file, hap1, cis2_file, hap2, trans_file, out_file):
    """
    Merges cis and trans pairs files, adding haplotype prefixes to specific columns.

    Args:
        cis1_file: Path to the first cis pairs file (gzipped).
        hap1: Haplotype prefix for the first cis pairs.
        cis2_file: Path to the second cis pairs file (gzipped).
        hap2: Haplotype prefix for the second cis pairs.
        trans_file: Path to the trans pairs file (gzipped).
        out_file: Path to the output merged pairs file.
    """
    try:
        with gzip.open(cis1_file, 'rt') as cis1, \
             gzip.open(cis2_file, 'rt') as cis2, \
             gzip.open(trans_file, 'rt') as trans, \
             open(out_file, 'w') as out:

            print(f"parsing {cis1_file}", file=sys.stderr)
            for line in cis1:
                if line.startswith('#'):
                    out.write(line)
                    continue
                a = line.strip().split('\t')
                a[1] += hap1
                a[3] += hap1
                out.write('\t'.join(a) + '\n')

            print(f"parsing {cis2_file}", file=sys.stderr)
            for line in cis2:
                if line.startswith('#'):
                    continue
                a = line.strip().split('\t')
                a[1] += hap2
                a[3] += hap2
                out.write('\t'.join(a) + '\n')

            print(f"parsing {trans_file}", file=sys.stderr)
            for line in trans:
                if line.startswith('#'):
                    continue
                a = line.strip().split('\t')
                a[1] += hap1
                a[3] += hap2
                out.write('\t'.join(a) + '\n')
        print("all done", file=sys.stderr)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("usage: merge_pairs.py cis_pairs1 prefix1 cis_pairs2 prefix2 trans_pairs out_pairs", file=sys.stderr)
        sys.exit(1)
    cis1_file, hap1, cis2_file, hap2, trans_file, out_file = sys.argv[1:]
    merge_pairs(cis1_file, hap1, cis2_file, hap2, trans_file, out_file)
