"""
The overwhelming majority of periodogram peaks are not real periods,
particularly near aliases of the TESS orbit.  In this script, we try to select
real astrophysical periods.

For every sector, do the following:

First, make a 2d histogram of the LSP peak value, bestlspval, versus
the period.  Save this as a plot, where the sector number is part of the plot
name.  Trim out any regions with N_COUNTS_CUTOFF of >= 5.

Next, using the output from the histogram calculation, plot the 1d distribution
of counts across all the bins.  Save this as a separate plot, similarly using
the sector number in the name.

Finally, write "0to500pc_quicklook_results_trimmed_periodcleaned" csv files,
including the "no duplicate, supplement" version.
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join
from trp.paths import RESULTSDIR
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap

def create_histograms(df, period_bin_edges, lsp_bin_edges, output_dir,
                      N_COUNTS_CUTOFF=5):
    """
    Create 2D and 1D histograms for each sector and save the plots.

    Args:
        df (pandas.DataFrame): DataFrame containing the data.
        period_bin_edges (numpy.ndarray): Bin edges for the period axis.
        lsp_bin_edges (numpy.ndarray): Bin edges for the LSP peak value axis.
        output_dir (str): Directory to save the output plots.
    """
    sectors = df['sector'].unique()
    selected_rows = {}

    for sector in np.sort(sectors):
        print(f"{sector}...")
        sector_df = df[df['sector'] == sector]

        # Create a 2D histogram of LSP peak value vs. period
        hist, _, _ = np.histogram2d(sector_df['period'], sector_df['bestlspval'],
                                    bins=[period_bin_edges, lsp_bin_edges])

        # Plot the 2D histogram
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.pcolormesh(period_bin_edges, lsp_bin_edges, hist.T,
                           vmin=1, vmax=N_COUNTS_CUTOFF, cmap='Blues')
        # Overlay gray "X" hatch on portions of the histogram with more than N counts
        masked_hist = np.ma.masked_where(hist.T <= N_COUNTS_CUTOFF, hist.T)
        cmap_gray = ListedColormap(['gray'])
        ax.pcolormesh(period_bin_edges, lsp_bin_edges, masked_hist,
                      cmap=cmap_gray, hatch='xxx', alpha=1)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Period')
        ax.set_ylabel('LSP Peak Value')
        ax.set_title(f'Sector {sector}: 2D Histogram of LSP Peak Value vs. Period')
        fig.colorbar(im, ax=ax, label='Counts')
        fig.savefig(join(output_dir, f'sector_{str(sector).zfill(4)}_2d_hist.png'))
        plt.close(fig)

        # Plot the 1D distribution of counts across all bins
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(hist.flatten(), bins=np.arange(0,1000+1,1))
        ax.set_xlabel('Counts')
        ax.set_ylabel('Frequency')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f'Sector {sector}: 1D Distribution of Counts')
        fig.savefig(join(output_dir, f'sector_{str(sector).zfill(4)}_1d_dist.png'))
        plt.close(fig)

        # Select rows with at most 2 counts in the histogram
        period_indices = np.digitize(sector_df['period'], period_bin_edges) - 1
        lsp_indices = np.digitize(sector_df['bestlspval'], lsp_bin_edges) - 1

        # Clip the indices to be within the valid range
        period_indices = np.clip(period_indices, 0, len(period_bin_edges) - 2)
        lsp_indices = np.clip(lsp_indices, 0, len(lsp_bin_edges) - 2)

        mask = hist[period_indices, lsp_indices] <= N_COUNTS_CUTOFF
        selected_rows[sector] = sector_df[mask]

    return selected_rows


def main():
    localdir = "/Users/luke/local/trp_RESULTS"
    # NOTE: includes duplicates, but that is OK for an exercise focused on removing
    # periods dominated by systematics.
    csvpath = join(localdir, "0to500pc_quicklook_results_trimmed.csv")
    output_dir = join(RESULTSDIR, 'systematics_cleaning')
    if not os.path.exists(output_dir): os.mkdir(output_dir)

    df = pd.read_csv(csvpath)
    df = df[~pd.isnull(df.period)]
    df['sector'] = df.sector.astype(int)

    period_bin_edges = np.logspace(np.log10(0.85), np.log10(13.5), 200)
    lsp_bin_edges = np.logspace(np.log10(1e-3), np.log10(1), 100)

    selected_rows = create_histograms(df, period_bin_edges, lsp_bin_edges, output_dir)

    # Print the number of selected rows for each sector
    N = 0
    sel_dfs = []
    for sector, sel_sector_df in selected_rows.items():
        print(f"Sector {sector}: {len(sel_sector_df)} rows selected")
        N += len(sel_sector_df)
        sel_dfs.append(sel_sector_df)
    print(f"Total: {N} (w/out LSPval cut)")

    selected_df = pd.concat(sel_dfs, ignore_index=True)
    sel = (
        (selected_df.bestlspval > 0.05)
        &
        (selected_df.period < 12)
    )
    selected_df = selected_df[sel]

    outcsvpath = join(
        localdir, "0to500pc_quicklook_results_trimmed_periodcleaned.csv"
    )
    selected_df.to_csv(
        outcsvpath, index=False
    )
    print(f'Wrote {outcsvpath}')

    from copy import deepcopy
    df = deepcopy(selected_df)
    if 'log10_ampl' not in df.columns:
        df['log10_ampl'] = np.log10(df['a_90_10_model'])
    if 'log10_period' not in df.columns:
        df['log10_period'] = np.log10(df['period'])

    gcsvpath = '/Users/luke/local/QLP/QLP_s1s55_X_GDR2_parallax_gt_2.csv'
    gdf = pd.read_csv(gcsvpath)
    selcols = 'dr2_source_id,ra,dec,parallax,pmra,pmdec,M_G,phot_g_mean_mag,bp_rp,ticid'.split(",")
    gdf = gdf[selcols]

    mdf = df.merge(gdf, on='ticid', how='inner')

    from sunnyhills.physicalpositions import calculate_XYZ_given_RADECPLX
    x,y,z = calculate_XYZ_given_RADECPLX(mdf.ra, mdf.dec, mdf.parallax)
    x0 = -8122
    mdf['x'] = x-x0
    mdf['y'] = y
    mdf['z'] = z

    suppcsvpath = join(localdir,
                       "0to500pc_quicklook_results_trimmed_periodcleaned_supp.csv")
    mdf.to_csv(suppcsvpath, index=False)
    print(f"Wrote {suppcsvpath}")

    mdf = mdf.sort_values(by='bestlspval', ascending=True)
    mdf = mdf.drop_duplicates(subset='ticid', keep='first')

    suppcsvpath = join(localdir,
                       "0to500pc_quicklook_results_trimmed_periodcleaned_supp_nodup.csv")
    mdf.to_csv(suppcsvpath, index=False)
    print(f"Wrote {suppcsvpath}")
    print(f"Total: {len(mdf)} (w/ LSPval cut, and P<12d, and no dups)")


if __name__ == '__main__':
    main()
