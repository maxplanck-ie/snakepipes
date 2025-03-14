#!/usr/bin/env python

import argparse
import logging
import pandas as pd

import cooler
import cooltools
import cooltools.lib.plotting
from cooltools import insulation

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.4'):
    raise AssertionError("tutorials rely on cooltools version 0.5.4 or higher,"+
                         "please check your cooltools version and update to the latest")



if __name__ == "__main__":
    # Defaults
    resolution = 10000
    threads = 1

    parser = argparse.ArgumentParser(description="Script to plot P(s) curves from mCool files")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-r", "--resolution", type=int, help=f"Base resolution, default={resolution}", default=resolution)
    parser.add_argument("-m", "--mcool", type=str, help="Path to mCool file", required=True)
    parser.add_argument("-o", "--output", type=str, help="Path to output file", required=True)
    parser.add_argument("-s", "--chromsizes", type=str, help="Path to chromomal sizes file", required=True)
    parser.add_argument("-c", "--chroms", type=str, help="Comma separated list of chromosomes", required=True)
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
        logging.debug("Verbose mode on")
    else:
        logging.basicConfig(level=logging.INFO)

    # get list of chromosomes
    chroms = args.chroms.split(',')
    logging.debug(f"Chromsomes: {chroms}")

    # load chromosomal sizes
    chromsizes = pd.read_csv(args.chromsizes,
                             sep='\t',
                             names=['chrom', 'size'],
                             dtype={'chrom': str, 'size': int})
    # filter chromosomal sizes
    chromsizes = chromsizes[chromsizes.chrom.isin(chroms)]
    logging.debug(f"Chromsizes: {chromsizes}")

    # define view area
    view_df = pd.DataFrame({'chrom': chromsizes['chrom'],
                            'start': 0,
                            'end': chromsizes['size'],
                            'name': chromsizes['chrom']})


    # load data
    logging.debug(f"Loading data from {args.mcool} with resolution {args.resolution}...")
    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{args.resolution}")
    
    resolution = args.resolution
    windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
    logging.debug(f"Windows: {windows}")

    logging.debug("computing insulation scores")
    insulation_table = insulation(clr, windows, view_df=view_df, verbose=True)

    logging.debug("writing final table")
    insulation_table.to_csv(args.output, sep = "\t")
