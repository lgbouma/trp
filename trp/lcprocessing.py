"""
Contents:
    | rotation_periodsearch
"""

#######################################
# ASTROBASE IMPORTS CAN BREAK LOGGING #
#######################################
from astrobase.periodbase import pgen_lsp, stellingwerf_pdm
from astrobase.checkplot import checkplot_png
from astrobase.lcmath import phase_magseries, phase_bin_magseries

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
import numpy as np, pandas as pd
from numpy import array as nparr

import os, multiprocessing, pickle

nworkers = multiprocessing.cpu_count()

def rotation_periodsearch(times, fluxs, starid, outdir, t0=None,
                          periodogram_method="ls", runperiodsearch=1,
                          do_finetune=0, write_pngs=0, nworkers=nworkers,
                          cachedict=None):
    """
    Given time and flux, calculate a periodogram.

    Plots and pickle files will be written to `outdir` using the `starid`
    string.

    Args:

        times (np.ndarray):
            Array of times.

        fluxs (np.ndarray):
            Array of fluxes.

        starid (str):
            Identifier used for cacheing.

        outdir (str):
            Path used for cacheing.

        t0 (None, str, int, or float):
            Epoch at which to phase.  None defaults to t0=1618.  Giving the
            string "binmin" defaults to phase-folding, and taking the
            arg-minimum.  Any int or float will be passed as the manual phase.

        periodogram_method (str):
            "pdm" (phase dispersion minimization) or "ls" (lomb-scargle).

        cachedict (dict or None):
            If passed, the keys and values in this dictionary will be passed to
            the pickle cache file.  This is relevant if for instance
            ```
                cachedict = {
                    'x_obs': x_obs,
                    'y_obs': y_obs
                }
            ```
            where x_obs and y_obs are some non-processed form of times and
            fluxs that you also may wish to visualize in the future.

    Returns:

        dict : results

            A dictionary of the results, containing:
                'lsp':periodogram results, 'fine_lsp':fine periodogram results,
                'times':times, 'fluxs':fluxs, 'period':fine_lsp['bestperiod'],
                't0': t0, 'outdir':outdir, 'periodogram_method': ...
            Note that just because the keys are "lsp", the actual method being
            used depends on periodogram_method
    """

    assert isinstance(starid, str)

    pklpath = os.path.join(outdir, f"{starid}_rotation_periodsearch.pkl")
    if os.path.exists(pklpath):
        LOGINFO(f"Found {pklpath}, loading and continuing.")
        with open(pklpath, 'rb') as f:
            d = pickle.load(f)
        return d

    sep = 1
    # eg., a single TESS sector at 2-minute cadence has 2e4 points.  this cuts
    # it down for period-search purposes to 4e3, which helps the runtime!
    if len(times) > 1e4:
        sep = 5
    # eg., a single TESS sector at 2-minute cadence has 1.2e5 points.
    if len(times) > 1e5:
        sep = 50

    startp, endp = 0.1, 10

    # for the fine-tuning
    delta_P = 0.2
    stepsize = 1e-5

    LOGINFO(f'Beginning period search for {starid}')

    if periodogram_method == 'ls' and runperiodsearch:

        lsp = pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0,
            nbestpeaks=10, nworkers=nworkers
        )

        if do_finetune:
            fine_lsp = pgen_lsp(
                times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
                startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
                endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
                autofreq=False, sigclip=5.0, stepsize=stepsize,
                nworkers=nworkers
            )
        else:
            LOGINFO(
                f"Found P={lsp['bestperiod']:.3f} d; skipping fine period "
                f"search."
            )
            fine_lsp = None

    elif periodogram_method == 'pdm' and runperiodsearch:

        lsp = stellingwerf_pdm(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0, nbestpeaks=10,
            nworkers=nworkers
        )

        if do_finetune:
            if lsp['bestperiod'] < 2:
                fine_lsp = stellingwerf_pdm(
                    times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
                    startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
                    endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
                    autofreq=False, sigclip=5.0, stepsize=stepsize,
                    nworkers=nworkers
                )
            else:
                LOGINFO(
                    f"Found P={lsp['bestperiod']:.3f} d; skipping fine period "
                    f"search."
                )
                fine_lsp = deepcopy(lsp)

    LOGINFO(42*'.')
    LOGINFO(f"Standard autofreq period: {lsp['bestperiod']:.7f} d")
    if do_finetune:
        LOGINFO(f"Fine period: {fine_lsp['bestperiod']:.7f} d")
        LOGINFO(f"Fine - standard: {fine_lsp['bestperiod']-lsp['bestperiod']:.7f} d")
    LOGINFO(42*'.')

    outfile = os.path.join(
        outdir, f'{starid}_{periodogram_method}_subset_checkplot.png'
    )
    if not os.path.exists(outfile) and write_pngs:
        try:
            checkplot_png(lsp, times, fluxs, fluxs*1e-4, magsarefluxes=True,
                          phasewrap=True, phasesort=True, phasebin=0.002,
                          minbinelems=7, plotxlim=(-0.6,0.6), plotdpi=75,
                          outfile=outfile, verbose=True)
        except Exception as e:
            LOGEXCEPTION(e)
            LOGINFO("Continuing...")
            pass

    if t0 is None:
        # default phase
        t0 = 1642.

    elif t0 == 'binmin':
        # bin the phase-fold to 50 points, take the minimum index.

        if fine_lsp is not None:
            period = fine_lsp['bestperiod']
        else:
            period = lsp['bestperiod']
        x,y = times, fluxs-np.nanmean(fluxs)
        t0_ini = np.nanmin(x)
        _pd = phase_magseries(x, y, period, t0_ini, wrap=False,
                              sort=False)
        x_fold = _pd['phase']
        y = _pd['mags']
        bs_days = period/50
        try:
            orb_bd = phase_bin_magseries(x_fold, y, binsize=bs_days, minbinelems=3)
            min_phase = orb_bd['binnedphases'][np.argmin(orb_bd['binnedmags'])]
            t0 = t0_ini + min_phase*period
        except (ValueError, TypeError):
            # can be raised for very short periods...
            t0 = t0_ini

    elif isinstance(t0, (int, float)):
        pass

    else:
        raise NotImplementedError

    if fine_lsp is not None:
        bestperiod = fine_lsp['bestperiod']
    else:
        bestperiod = lsp['bestperiod']
    d = {
        'lsp':lsp,
        'fine_lsp':fine_lsp,
        'times':times,
        'fluxs':fluxs,
        'period':bestperiod,
        't0': t0,
        'outdir':outdir,
        'periodogram_method': periodogram_method
    }
    if cachedict is not None:
        for k,v in cachedict.items():
            d[k] = v

    with open(pklpath, 'wb') as f:
        pickle.dump(d, f)
        LOGINFO(f'Made {pklpath}')

    return d
