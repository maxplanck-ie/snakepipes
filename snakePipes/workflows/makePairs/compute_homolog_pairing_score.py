#!/usr/bin/env python

import warnings
import numpy as np
import scipy.stats as st
import pandas as pd
import cooler
import cooltools
import cooltools.lib
from cooltools.lib import numutils

import click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('COOLER_PATH', nargs=1)
@click.option(
    '-o',
    '--output',
    type=click.File('w'),
    default='-',
    show_default=True,
    help='path to the output bedgraph file (by default, printed into stdout)')
@click.option(
    '-r',
    '--resolution',
    type=int,
    default=10000,
    show_default=True,
    help='base pair resolution of the cooler file')
@click.option(
    '--window-bp',
    type=int,
    default=20000, 
    show_default=True,
    help='window size to aggregate the contact frequency between homologous loci')
@click.option(
    '--homolog-suffixes', 
    nargs=2,
    default = ('_057', '_439'),
    show_default=True,
    help='suffixes used to distinguish the pair of homolog chromosomes')
@click.option(
    '--ignore-diags',
    type=int,
    default=0,
    show_default=True,
    help='the number of diagonals of the Hi-C matrix to ignore due to '
         'contamination by short-distance Hi-C ')
@click.option(
    '--poisson-perc',
    type=float,
    default=None,
    show_default=True,
    help='if provided, report the N-th percentile of the sampling distribution '
        'of the pairing score at each location, assuming the trans-homolog '
        'interaction counts are Poisson-distributed')
@click.option(
    '--balance',
    type=bool,
    default=True,
    show_default=True,
    help='use interaction frequencies normalized by iterative correction')
@click.option(
    '--normalize-by-cis',
    type=click.Choice(['False', 'median', 'True',]),
    default='False',
    show_default=True,
    help='divide the pairing score by the cis interaction frequency in the corresponding window. '
    'If median, divide the pairing score by the median of the cis interaction frequencies.')
@click.option(
    '--normalize-by-median',
    type=bool,
    default=True,
    show_default=True,
    help='divide the pairing score by its genome-wide median')
@click.option(
    '--transform',
    type=click.Choice(['linear', 'log2', 'log10']),
    default='log2',
    show_default=True,
    help='report either log-transformed or non-transformed (linear) pairing scores')
@click.option(
    '--report-per-homolog',
    type=bool,
    default=False,
    show_default=True,
    help='if True, duplicate the pairing scores for each of the homologs '
         'otherwise, report the pairing score once for each homolog pair'
    )
@click.option(
    '--agg-func',
    type=click.Choice(['nansum', 'nanmean', 'nanmedian']),
    default='nanmean',
    show_default=True,
    help='the function to calculate score over diagonal pixels')
@click.option(
    '--as-bedgraph',
    is_flag=True,
    show_default=True,
    help='output a bedgraph with the pairing score instead of a full table')

def make_pairing_bedgraph(
    cooler_path,
    resolution,
    output,
    window_bp,
    homolog_suffixes,
    ignore_diags,
    poisson_perc,
    balance,
    normalize_by_cis,
    normalize_by_median,
    transform,
    report_per_homolog, 
    agg_func,
    as_bedgraph
):
    """Generate a bedgraph file with a homolog pairing score.
    """

    clr = cooler.Cooler(f"{cooler_path}::/resolutions/{resolution}")
    
    df = get_homolog_pairing_score(
        clr,
        window_bp,
        homolog_suffixes,
        ignore_diags,
        poisson_perc,
        balance,
        normalize_by_cis,
        report_per_homolog,
        agg_func
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if transform == 'log2':
            df['pairing'] = np.log2(df['pairing'])
            df['pairing'].loc[~np.isfinite(df['pairing'])] = np.nan
            df['cis_norm'] = np.log2(df['cis_norm'])
            df['cis_norm'].loc[~np.isfinite(df['cis_norm'])] = np.nan
        elif transform == 'log10':
            df['pairing'] = np.log10(df['pairing'])
            df['pairing'].loc[~np.isfinite(df['pairing'])] = np.nan
            df['cis_norm'] = np.log10(df['cis_norm'])
            df['cis_norm'].loc[~np.isfinite(df['cis_norm'])] = np.nan

    gw_median_pairing = np.nanmedian(df['pairing'])
    gw_median_cis_norm = np.nanmedian(df['cis_norm'])
    if normalize_by_median:
        if transform.startswith('log'): 
            df['pairing'] -= gw_median_pairing
            df['cis_norm'] -= gw_median_cis_norm
        else:
            df['pairing'] /= gw_median_pairing
            df['cis_norm'] /= gw_median_cis_norm

    if as_bedgraph:
        df['pairing'] = np.nan_to_num(df['pairing'])
        df['cis_norm'] = np.nan_to_num(df['cis_norm'])
        
    bedgraph_header = (
        """track type=bedGraph
        name="{transform}pairing{percentile}"
        description="{transform}pairing{percentile}, {window_bp}bp window at {res}bp resolution{norm_by_cis}{normalize_by_median}{ignore_diags}. 
        Summary: trans median: {trans_median}, cis median: {cis_median}"
        visibility=full color=0,0,0 altColor=50,50,50 graphType=points yLineMark=0 yLineOnOff=on
        group="{transform}pairing"
        """.replace('\n', ' ').format(
            transform = (transform +'_') if transform else '',
            percentile = ('_'+str(poisson_perc)+'th_perc') if poisson_perc else '',
            window_bp=window_bp,
            res = clr.info['bin-size'],
            norm_by_cis = {'True':', normalized by cis diagonal',
                           'median':', normalized by cis diagonal median',
                           'False':'',
                            }[normalize_by_cis],
            normalize_by_median = ', normalized by the GW-median' if normalize_by_median else '',
            ignore_diags = ', ignored {} diagonals'.format(ignore_diags) if ignore_diags else '',
            trans_median = gw_median_pairing,
            cis_median = gw_median_cis_norm
            )
        + '\n')
    

    if as_bedgraph:
        output.write(bedgraph_header)

        df.to_csv(
            output.iloc[:, :4],
            sep='\t',
            header=None,
            index=False)
    else:
        output.write('#' + bedgraph_header)

        df.to_csv(
            output,
            sep='\t',
            index=False)


def get_homolog_pairing_score(
    clr,
    window_bp=10000,
    homolog_suffixes = ('_057', '_439'),
    ignore_diags = 0,
    poisson_perc = None,
    balance = True,
    normalize_by_cis = 'False',
    report_per_homolog = True,
    agg_func='nanmean'
):
    """Compute the homolog pairing score
    """
    pairing_dfs = []

    hom_chroms = clr.chromnames
    base_chroms = []
    if homolog_suffixes == ('',''):
        base_chroms = hom_chroms
    else:
        for chrom in hom_chroms:
            base_chrom = ''
            
            if chrom.endswith(homolog_suffixes[0]):
                base_chrom = chrom[:-len(homolog_suffixes[0])]
            if chrom.endswith(homolog_suffixes[1]):
                base_chrom = chrom[:-len(homolog_suffixes[1])]
            if (
                ((base_chrom+homolog_suffixes[0]) in hom_chroms)
                and ((base_chrom+homolog_suffixes[1]) in hom_chroms)
                and base_chrom not in base_chroms
            ):

                base_chroms.append(base_chrom)

    window = int(np.floor(window_bp / 2 / clr.info['bin-size']))
    bins = clr.bins()[:]
    agg_func = {'nanmean':np.nanmean, 
                'nansum':np.nansum,
                'nanmedian':np.nanmedian}[agg_func]

    for chrom in base_chroms:
        trans_diag_raw = None
        trans_diag_balanced = None
        if (not balance) or (poisson_perc is not None):
            trans_mat_raw = clr.matrix(balance=False).fetch(
                chrom + homolog_suffixes[0], chrom+homolog_suffixes[1]).astype(float)
            trans_diag_raw = _take_big_diagonal_pixel(
                trans_mat_raw, window, ignore_diags, agg_func=agg_func)
        if balance:
            trans_mat_balanced = clr.matrix(balance=True).fetch(
                chrom + homolog_suffixes[0], chrom+homolog_suffixes[1]).astype(float)
            trans_diag_balanced = _take_big_diagonal_pixel(
                trans_mat_balanced, window, ignore_diags, agg_func=agg_func)
           
        trans_diag = trans_diag_balanced if balance else trans_diag_raw
        
        if poisson_perc is not None:
            trans_diag *= st.poisson.isf(poisson_perc/100, trans_diag_raw) / trans_diag_raw

        cis_hom1_mat = clr.matrix(balance=balance).fetch(
            chrom + homolog_suffixes[0], chrom+homolog_suffixes[0]).astype(float)
        cis_hom2_mat = clr.matrix(balance=balance).fetch(
            chrom + homolog_suffixes[1], chrom+homolog_suffixes[1]).astype(float)
        cis_hom1_diag = _take_big_diagonal_pixel(
            cis_hom1_mat, window, ignore_diags, agg_func=agg_func)
        cis_hom2_diag = _take_big_diagonal_pixel(
            cis_hom2_mat, window, ignore_diags, agg_func=agg_func)

        chrom_bins = bins[bins.chrom == chrom+homolog_suffixes[0]]
        pairing_df = pd.DataFrame(dict(binid=chrom_bins.index.values))
        pairing_df['chrom'] = chrom
        pairing_df['start'] = chrom_bins.start.values
        pairing_df['end'] = chrom_bins.end.values
        if normalize_by_cis == 'True':
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pairing_df['pairing'] = trans_diag / ((cis_hom1_diag + cis_hom2_diag)/2)
        elif normalize_by_cis == 'median':
            pairing_df['pairing'] = trans_diag / np.nanmedian((cis_hom1_diag + cis_hom2_diag)/2)
        elif normalize_by_cis == 'False':
            pairing_df['pairing'] = trans_diag
        else:
            raise ValueError(
                'normalize_by_cis can only take string values of "True", "False" and "median"')
        pairing_df['cis_norm1'] = cis_hom1_diag
        pairing_df['cis_norm2'] = cis_hom2_diag

        if report_per_homolog:
            pairing_df_1 = pairing_df.copy()
            pairing_df_2 = pairing_df.copy()
            pairing_df_1['chrom'] = chrom + homolog_suffixes[0]
            pairing_df_2['chrom'] = chrom + homolog_suffixes[1]
            pairing_df_1['binid'] = bins[
                bins['chrom'] == chrom + homolog_suffixes[0]].index.values
            pairing_df_2['binid'] = bins[
                bins['chrom'] == chrom + homolog_suffixes[1]].index.values
            pairing_df_1['cis_norm'] = cis_hom1_diag
            pairing_df_2['cis_norm'] = cis_hom2_diag

            if (np.searchsorted(bins['chrom'], chrom + homolog_suffixes[0], 'left') <
                np.searchsorted(bins['chrom'], chrom + homolog_suffixes[1], 'left')):

                pairing_dfs.append(pairing_df_1)
                pairing_dfs.append(pairing_df_2)
            else:
                pairing_dfs.append(pairing_df_2)
                pairing_dfs.append(pairing_df_1)
        else:
            pairing_dfs.append(pairing_df)


    pairing_dfs = pd.concat(pairing_dfs, ignore_index=True)
    pairing_dfs.sort_values('binid', inplace=True)
    pairing_dfs.drop('binid', axis=1, inplace=True)

    return pairing_dfs

def _take_big_diagonal_pixel(
    mat,
    pad_bins=10,
    ignore_diags=3,
    agg_func=np.nanmean
):
    """Operations over the main diagonal values
    """
    if (ignore_diags):
        mat = mat.copy()
        for i in range(-ignore_diags, ignore_diags+1):
            numutils.set_diag(mat, np.nan, i)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        N = mat.shape[0]
        score = np.nan * np.ones(N)
        for i in range(0, N):
            lo = max(0, i-pad_bins)
            hi = min(i+pad_bins+1, N)
            # nanmean of interactions to reduce the effect of bad bins
            score[i] = agg_func(mat[lo:hi, lo:hi])

    return score

if __name__ == '__main__':
    make_pairing_bedgraph()
