#!/usr/bin/env python

# import core packages
import logging
import argparse

# import additional packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool

# import open2c libraries
import cooler
import cooltools


if __name__ == "__main__":
    # Defaults
    resolution = 10000
    threads = 1

    parser = argparse.ArgumentParser(description="Script to plot P(s) curves from mCool files")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-r", "--resolution", type=int, help=f"Base resolution, default={resolution}", default=resolution)
    parser.add_argument("-m", "--mcool", type=str, help="Path to mCool file", required=True)
    parser.add_argument("-o", "--output", type=str, help="Path to output file", required=True)
    parser.add_argument("-t", "--threads", type=int, help=f"Number of threads, default={threads}", default=threads)
    parser.add_argument("-s", "--chromsizes", type=str, help="Path to chromosome sizes file", required=True)
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
                             dtype={'chrom': str, 'size': int}
                            )
    # filter chromosomal sizes
    chromsizes = chromsizes[chromsizes.chrom.isin(chroms)]
    logging.debug(f"Chromsizes: {chromsizes}")
    
    # load data
    logging.debug(f"Loading data from {args.mcool} with resolution {args.resolution}...")
    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{args.resolution}")

    # compute P(s) curves
    logging.debug("Computing P(s) curves...")
    cvd_smooth_agg = cooltools.expected_cis(
        clr=clr,
        view_df=None,
        smooth=True,
        aggregate_smoothed=True,
        smooth_sigma=0.1,
        nproc=threads
    )

    # filter short distances
    logging.debug("Filtering short distances...")
    cvd_smooth_agg['balanced.avg.smoothed'].loc[cvd_smooth_agg['dist'] < 2] = np.nan

    # plot
    logging.debug("Plotting P(s) curves...")
    fig, ax = plt.subplots(1,1)

    for region in chromsizes['chrom']:
        ax.loglog(
            cvd_smooth_agg['dist_bp'].loc[cvd_smooth_agg['region1']==region],
            cvd_smooth_agg['balanced.avg.smoothed'].loc[cvd_smooth_agg['region1']==region],
        )
        ax.set(
            xlabel='separation, bp',
            ylabel='IC contact frequency')
        ax.set_aspect(1.0)
        ax.grid(lw=0.5)

    ax.legend(chromsizes['chrom'], loc='upper right')
    fig.savefig(args.output)
    plt.close(fig)

    # saving data table
    logging.debug("Saving data table...")
    cvd_smooth_agg.to_csv(args.output.replace('.png', '.tsv'), sep='\t') 