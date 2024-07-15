"""
Contents:
    | plot_rotvetter
"""
#######################################
# ASTROBASE IMPORTS CAN BREAK LOGGING #
#######################################
from astrobase.lcmath import (
    phase_magseries, phase_bin_magseries, sigclip_magseries,
    find_lc_timegroups, phase_magseries_with_errs, time_bin_magseries
)
from astrobase.services.identifiers import tic_to_gaiadr2, tic_to_simbad
from astrobase.lcmath import sigclip_magseries

#############
## LOGGING ##
#############
import logging
from trp import log_sub, log_fmt, log_date_fmt

DEBUG = False
if DEBUG:
    level = logging.DEBUG
else:
    level = logging.INFO
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=level,
    style=log_sub,
    format=log_fmt,
    datefmt=log_date_fmt,
    force=True
)

LOGDEBUG = LOGGER.debug
LOGINFO = LOGGER.info
LOGWARNING = LOGGER.warning
LOGERROR = LOGGER.error
LOGEXCEPTION = LOGGER.exception

#############
## IMPORTS ##
#############
import os, pickle
from os.path import join
from glob import glob
from datetime import datetime
from copy import deepcopy
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from numpy import array as nparr
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams

from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

import matplotlib.patheffects as pe
from matplotlib.ticker import MaxNLocator, FixedLocator, FuncFormatter
from matplotlib.transforms import blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator
)


from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d, PchipInterpolator

from aesthetic.plot import savefig, format_ax, set_style
from astroquery.mast import Catalogs
from complexrotators.plotting import (
    plot_phased_light_curve, get_ylimguess, get_2px_neighbors, dss_overlay
)
from cdips.utils.gaiaqueries import (
    given_dr2_sourceids_get_edr3_xmatch, given_source_ids_get_gaia_data
)

def plot_rotvetter(
    outpath,
    lcpath,
    starid,
    periodsearch_result=None,
    binsize_phase=0.005,
    lcpipeline='qlp'
    ):
    """
    periodsearch_result contains:
    """

    # get data
    d = periodsearch_result
    d['times'], d['fluxs'] = d['x_obs'], d['y_obs']
    hdul = fits.open(lcpath)
    hdr = hdul[0].header
    data = hdul[1].data
    hdul.close()

    # quality flags
    QUALITYKEYDICT = {
        'spoc2min': 'QUALITY',
        'qlp': 'QUALITY',
        'cdips': 'IRQ3'
    }
    qual = data[QUALITYKEYDICT[lcpipeline]]
    if lcpipeline in ['spoc2min', 'qlp']:
        sel = (qual == 0)
    elif lcpipeline == 'cdips':
        sel = (qual == 'G')

    # centroid data (read; not used)
    CENTRKEYDICT = {
        'spoc2min': [
            'MOM_CENTR2', # column
            'MOM_CENTR1' # row
        ],
        'qlp': [
            'SAP_X', # column
            'SAP_Y' # row
        ],
        'cdips': [
            'XIC', # column
            'YIC' # row
        ]
    }
    xc = data[CENTRKEYDICT[lcpipeline][0]][sel] # column
    yc = data[CENTRKEYDICT[lcpipeline][1]][sel] # row

    # background data
    BKGDKEYDICT = {
        'spoc2min': 'SAP_BKG',
        'qlp': 'SAP_BKG',
        'cdips': 'BGV'
    }
    bgv = data[BKGDKEYDICT[lcpipeline]][sel]

    assert len(xc) == len(d['times'])

    # make plot
    plt.close('all')
    set_style("clean")

    #fig = plt.figure(figsize=(8,4.5))
    fig = plt.figure(figsize=(8,6))
    axd = fig.subplot_mosaic(
        """
        AAAAAACC
        AAAAAACC
        BBBDDDEE
        BBBDDDEE
        FFFFFFEE
        """
    )

    axd['A'].get_shared_x_axes().join(axd['A'], axd['F'])

    #
    # sap flux vs time (2hr bin over whatever cadence)
    #
    ax = axd['A']

    bd = time_bin_magseries(d['times'], d['fluxs'], binsize=7200, minbinelems=1)
    yoffset = np.nanmean(bd['binnedmags'])
    ax.scatter(d['times'], 1e2*(d['fluxs']-yoffset),
               c='lightgray', s=0.5, zorder=1)
    ax.scatter(bd['binnedtimes'], 1e2*(bd['binnedmags']-yoffset),
               c='k', s=3, zorder=2, rasterized=True)

    yval = np.nanpercentile(1e2*(bd['binnedmags']-yoffset), 95)

    _time, _flux = bd['binnedtimes'], 1e2*(bd['binnedmags']-yoffset)
    x0 = _time[np.argmin(np.abs(_flux - np.nanpercentile(_flux, 1)))]
    x1 = x0 + d['period']
    xval = x0 + (x1 - x0)/2
    xerr = x1 - xval
    #ax.hlines(yval, ymin, ymax, colors='red', alpha=1,
    #          linestyles='-', zorder=5, linewidths=2)
    ax.errorbar(xval, yval, xerr=xerr, alpha=1,
                marker='.', elinewidth=1, capsize=2, lw=0, mew=0.1,
                color='red', markersize=0, zorder=5)

    ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))
    ax.update({'xlabel': 'Time [BTJD]', 'ylabel': 'SAP Flux [%]', 'ylim': ylim})

    #
    # pdm periodogram
    #
    ax = axd['B']

    ax.plot(d['lsp']['periods'], d['lsp']['power'], c='k', lw=1)
    ax.scatter(d['lsp']['nbestperiods'][:5], d['lsp']['nbestlspvals'][:5],
               marker='v', s=5, linewidths=0, edgecolors='none',
               color='C0', alpha=0.5, zorder=1000)
    ymin, ymax = ax.get_ylim()
    ax.vlines(d['period'], ymin, ymax, colors='darkgray', alpha=1,
              linestyles='-', zorder=-2, linewidths=1)
    P_harmonics = []
    for ix in range(1,11):
        P_harmonics.append(ix*d['period'])
        P_harmonics.append(d['period']/ix)

    sel = (
        (P_harmonics > np.nanmin(d['lsp']['periods']))
        &
        (P_harmonics < np.nanmax(d['lsp']['periods']))
    )
    ax.vlines(nparr(P_harmonics)[sel], ymin, ymax, colors='darkgray',
              alpha=0.5, linestyles=':', zorder=-2, linewidths=0.5)
    ax.set_ylim([ymin, ymax])
    ax.update({'xlabel': 'Period [d]', 'ylabel': 'LS Power', 'xscale': 'log'})

    #
    # phased LC
    #
    ax = axd['D']
    ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
        xlim=[-0.6,0.6], binsize_phase=0.01, BINMS=7, titlestr=None,
        showtext=False, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.9, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=True
    )
    txt = f'{d["period"]:.1f} d'
    ax.text(
        0.95, 0.05, txt, transform=ax.transAxes, fontsize='medium', ha='right',
        va='bottom'
    )
    ax.set_ylabel("$\Delta$ Flux [%]")
    ax.set_xticklabels(['-0.5', '', '0', '', '0.5'])
    ax.set_xlabel("Phase, φ")

    #
    # flux in BGD aperture
    #
    ax = axd['F']
    nbgv = bgv/np.nanmedian(bgv)

    bd = time_bin_magseries(d['times'], nbgv, binsize=7200, minbinelems=1)
    yoffset = np.nanmean(bd['binnedmags'])
    ax.scatter(d['times'], nbgv, c='lightgray', s=0.5, zorder=1)
    ax.scatter(bd['binnedtimes'], bd['binnedmags'],
               c='k', s=3, zorder=2, rasterized=True)
    ylim = get_ylimguess(bd['binnedmags'])
    ax.update({'xlabel': 'Time [BTJD]', 'ylabel': 'BGV/med(BGV)', 'ylim': ylim})

    # join subplots

    #
    # # star info (Gaia, TIC8, dip search)
    #
    ax = axd['E']
    ax.set_axis_off()

    # tic8 info
    TEFFKEYDICT = {
        'spoc2min': 'TEFF',
        'qlp': 'TEFF',
        'cdips': 'TICTEFF'
    }
    PMRAKEYDICT = {
        'spoc2min': 'PMRA',
        'qlp': 'PMRA',
        'cdips': 'PM_RA[mas/yr]'
    }
    PMDECKEYDICT = {
        'spoc2min': 'PMDEC',
        'qlp': 'PMDEC',
        'cdips': 'PM_Dec[mas/year]'
    }

    ticid = str(hdr['TICID'])
    sector = str(hdr['SECTOR'])
    cam = str(hdr['CAMERA'])
    ccd = str(hdr['CCD'])
    Tmag = f"{hdr['TESSMAG']:.1f}"
    if hdr[TEFFKEYDICT[lcpipeline]] is not None and hdr[TEFFKEYDICT[lcpipeline]] != 'nan':
        teff_tic8 = f"{int(hdr[TEFFKEYDICT[lcpipeline]]):d} K"
    else:
        teff_tic8 = f"NaN"
    if hdr['RA_OBJ'] != 'nan':
        ra = f"{hdr['RA_OBJ']:.2f}"
    else:
        ra = 'NaN'
    if hdr['DEC_OBJ'] != 'nan':
        dec = f"{hdr['DEC_OBJ']:.2f}"
    else:
        dec = 'NaN'
    if hdr[PMRAKEYDICT[lcpipeline]] is not None and hdr[PMRAKEYDICT[lcpipeline]] != 'nan':
        pmra = f"{hdr[PMRAKEYDICT[lcpipeline]]:.1f}"
    else:
        pmra = 'NaN'
    if hdr[PMDECKEYDICT[lcpipeline]] is not None and hdr[PMDECKEYDICT[lcpipeline]] != 'nan':
        pmdec = f"{hdr[PMDECKEYDICT[lcpipeline]]:.1f}"
    else:
        pmdec = 'NaN'
    ra_obj, dec_obj = hdr['RA_OBJ'], hdr['DEC_OBJ']
    c_obj = SkyCoord(ra_obj, dec_obj, unit=(u.deg), frame='icrs')

    # gaia info
    dr2_source_id = tic_to_gaiadr2(ticid)
    simbad_name = tic_to_simbad(ticid)

    runid = f"dr2_{dr2_source_id}"

    dr2_source_ids = np.array([np.int64(dr2_source_id)])
    try:
        gdf = given_source_ids_get_gaia_data(
            dr2_source_ids, runid, n_max=5, overwrite=False,
            enforce_all_sourceids_viable=True, savstr='', which_columns='*',
            table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=False
        )
    except Exception as e:
        LOGERROR(f'{ticid} = {runid} failed due to:\n{e}...')
        return 1
    try:
        gdf_ruwe = given_source_ids_get_gaia_data(
            dr2_source_ids, runid+"_ruwe", n_max=5, overwrite=False,
            enforce_all_sourceids_viable=True, savstr='', which_columns='*',
            table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=True
        )
        ruwe =  f"{gdf_ruwe['ruwe'].iloc[0]:.2f}"
    except AttributeError:
        ruwe = 'NaN'

    Gmag = f"{gdf['phot_g_mean_mag'].iloc[0]:.1f}"
    Rpmag = f"{gdf['phot_rp_mean_mag'].iloc[0]:.1f}"
    Bpmag = f"{gdf['phot_bp_mean_mag'].iloc[0]:.1f}"
    bp_rp = f"{gdf['bp_rp'].iloc[0]:.2f}"
    plx = f"{gdf['parallax'].iloc[0]:.2f}"
    plx_err = f"{gdf['parallax_error'].iloc[0]:.2f}"
    dist_pc = 1/(float(plx)*1e-3)
    dist = f"{dist_pc:.1f}"

    # nbhr info
    ticids, tmags = get_2px_neighbors(c_obj, hdr['TESSMAG'])
    brightest_inds_first = np.argsort(tmags)
    ticids = ticids[brightest_inds_first]
    tmags = tmags[brightest_inds_first]
    N_nbhrs = len(ticids)-1
    nbhrstr = ''
    MAX_N = 7
    if N_nbhrs >= 1:
        for _ticid, tmag in zip(ticids[1:MAX_N], tmags[1:MAX_N]):
            nbhrstr += f"TIC {_ticid}: ΔT={tmag-hdr['TESSMAG']:.1f}\n"

    txt = (
        f"{simbad_name}\n"
        f"TIC {ticid}\n"
        f"GDR2 {dr2_source_id}\n"
        f"SEC{sector}, CAM{cam}, CCD{ccd}\n"
        "—\n"
        f"α={ra}, δ={dec} (deg)\n"
        f"T={Tmag}\n"
        f"TIC8 Teff={teff_tic8}\n"
        f"G={Gmag}, RP={Rpmag}, BP={Bpmag}\n"
        f"BP-RP={bp_rp}\n"
        f"RUWE={ruwe}\n"
        f"plx={plx}"+"$\pm$"+f"{plx_err} mas\n"
        f"d={dist} pc\n"
        "—\n"
        f'A={1e2*d["a_90_10_model"]:.1f}%\n'
        f'P2P={1e2*d["p2p_noise"]:.1f}%\n'
        f'SN={d["snr_metric"]:.1f}\n'
        "—\n"
        'N$_{\mathrm{nbhr}}$: '+f'{N_nbhrs}\n'
        f'{nbhrstr}\n'
    )

    txt_x = -0.1
    txt_y = 0.5
    ax.text(txt_x, txt_y, txt, ha='left', va='center', fontsize='small', zorder=2,
            transform=ax.transAxes)

    #
    # DSS query
    #
    ra = hdr['RA_OBJ']
    dec = hdr['DEC_OBJ']

    dss_overlay(fig, axd, ra, dec)

    # set naming options
    s = ''

    # height / weight
    fig.tight_layout(h_pad=0.3)

    savefig(fig, outpath)
