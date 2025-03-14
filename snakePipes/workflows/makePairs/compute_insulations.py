#!/usr/bin/env python

# import core packages
import logging
import argparse

# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import python package for working with cooler files and tools for analysis
import cooler
from cooltools import insulation

if __name__ == "__main__":
    # Defaults
    resolution = 10000
    windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
    parser = argparse.ArgumentParser(description="Script to plot P(s) curves from mCool files")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-r", "--resolution", type=int, help=f"Base resolution, default={resolution}", default=resolution)
    parser.add_argument("-m", "--mcool", type=str, help="Path to mCool file", required=True)
    parser.add_argument("-o", "--output", type=str, help="Path to output file", required=True)
    parser.add_argument("-w", "--windows", type=str, help="Comma separated list of resolution windows, default={windows}", default=windows)
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
        logging.debug("Verbose mode on")
    else:
        logging.basicConfig(level=logging.INFO)

    # load data
    logging.debug(f"Loading data from {args.mcool} with resolution {args.resolution}...")
    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{args.resolution}")

    # compute insulations
    logging.debug("Computing insulations...")
    insulation_table = insulation(clr, windows, verbose=True)

    for window in windows:
        logging.debug(f"Writing {window} insulations to {args.output + f'_w{window}.tsv'}")
        window_summary =insulation_table.columns[[ str(window) in i for i in insulation_table.columns]]

        ins_table = insulation_table[['chrom','start','end','region','is_bad_bin'] + list(window_summary)].iloc[1000:1005]
        ins_table.to_csv(args.output + f'_w{window}.tsv', sep='\t', index=False)

    
    histkwargs = dict(
        bins=10**np.linspace(-4,1,200),
        histtype='step',
        lw=2,
    )

    fig, axs = plt.subplots(len(windows),1, sharex=True, figsize=(6,6), constrained_layout=True)
    for i, (w, ax) in enumerate(zip(windows, axs)):
        ax.hist(
            insulation_table[f'boundary_strength_{w}'],
            **histkwargs
        )
        ax.text(0.02, 0.9,
                f'Window {w//1000}kb',
                ha='left',
                va='top',
                transform=ax.transAxes)

        ax.set(
            xscale='log',
            ylabel='# boundaries'
        )

    axs[-1].set(xlabel='Boundary strength')
    fig.savefig(args.output + '.png', dpi=300)
    plt.close(fig)