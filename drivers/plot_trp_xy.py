import os
from os.path import join
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from numpy import array as nparr
import matplotlib as mpl

from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

from aesthetic.plot import savefig, format_ax, set_style

def plot_trp_xy(noaxis=1, dcut=500, source='bouma_trp'):

    set_style('clean')

    if source == 'bouma_trp':
        csvkey = '0to500pc_quicklook_results_trimmed_periodcleaned_supp_nodup_dr3.csv'
        csvpath = join('/Users/luke/local/trp_RESULTS', csvkey)
    elif source == 'boyle_tars':
        #csvkey = 'all_rotations_with_gaia_data_oct14_powerGTp1_Plt10_utoppower_supp.csv'
        csvkey = 'all_rotations_with_gaia_data_oct14_powerGTp1_periodLT10_some_period_cleaning.csv'
        csvpath = join('/Users/luke/Dropbox/proj/kairos_plato/data', csvkey)
    df = pd.read_csv(csvpath)

    x,y,z = nparr(df.x), nparr(df.y), nparr(df.z)

    x_sun, y_sun = 0, 0
    assert isinstance(dcut, int)

    r = np.sqrt( (x-x_sun)**2 + y**2 + z**2 )

    sel = (
        (r < dcut)
        &
        (df.ruwe < 1.4)
        #&
        #(~df.non_single_star)
    )

    x,y,z = x[sel], y[sel], z[sel]

    cutevery = 1
    x,y,z = x[::cutevery], y[::cutevery], z[::cutevery]

    fig, ax = plt.subplots(figsize=(4,4))

    showsun = 0
    if showsun:
        ax.scatter(
            x_sun, y_sun, c='black', alpha=1, zorder=1, s=20, rasterized=True,
            linewidths=1, marker='x'
        )

    # By default, just show all the stars as the same color.  The
    # "rasterized=True" kwarg here is good if you save the plots as pdfs,
    # to not need to save the positions of too many points.
    if dcut == 500:
        s = 0.3
        alpha = 0.3
    elif dcut == 250:
        s = 1
        alpha = 0.2

    ax.scatter(
        x, y, c='black', alpha=alpha, zorder=2, s=s, rasterized=True,
        linewidths=0, marker='.'
    )

    if noaxis:
        ax.set_axis_off()
    else:
        ax.set_xlabel("X [pc]")
        ax.set_ylabel("Y [pc]")

    s = f'_dcut{dcut}'
    if noaxis:
        s += f"_noaxis"

    outdir = '../results/'+csvkey.rstrip(".csv")
    if not os.path.exists(outdir): os.mkdir(outdir)

    outpath = join(outdir, f'trp_XY{source}_{s}.png')
    fig.savefig(outpath, bbox_inches='tight', dpi=400)
    print(f"Made {outpath}")

if __name__ == "__main__":
    plot_trp_xy(noaxis=1, dcut=250, source='boyle_tars')
    plot_trp_xy(noaxis=1, dcut=500, source='boyle_tars')
    assert 0
    plot_trp_xy(noaxis=1, dcut=250)
    plot_trp_xy(noaxis=1, dcut=500)
