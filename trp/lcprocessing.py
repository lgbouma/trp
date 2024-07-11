"""
Contents:
    | rotation_periodsearch
    | calculate_lsp
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

from astropy.timeseries import LombScargle

import os, multiprocessing, pickle

from trp.lcutils import p2p_rms

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


def calculate_lsp(times, fluxes, starid=None, outdir=None, nyquist_factor=2,
                  samples_per_peak=10,
                  n_best=5, period_min=0.1, period_max=27/2, N_freq=int(1e6),
                  cachedict=None, periodogram_method='astropyls',
                  periodepsilon=0.1, cache_periodogram_pkls=1):
    """Calculate the Lomb Scargle periodogram for the given times and fluxes.

    Args:
        times (numpy.ndarray): Array of time values.
        fluxes (numpy.ndarray): Array of flux values.
        starid (str): Identifier used for cacheing.
        outdir (str): Path used for cacheing.
        nyquist_factor (float, optional): Oversampling factor for the frequency
            grid. Defaults to 2.
        n_best (int, optional): Number of best periods to return. Defaults to 5.

    Returns:
        dict: A dictionary containing the periodogram results, peak period, and
            top N best periods and values.
    """

    assert isinstance(starid, str)

    pklpath = os.path.join(outdir, f"{starid}_rotation_periodsearch.pkl")
    if os.path.exists(pklpath):
        LOGINFO(f"Found {pklpath}, loading and continuing.")
        with open(pklpath, 'rb') as f:
            d = pickle.load(f)
        return d

    # Create a LombScargle object
    ls = LombScargle(times, fluxes)

    # Calculate the frequency grid
    minimum_frequency = 1/period_max
    maximum_frequency = 1/period_min
    frequency = np.linspace(minimum_frequency, maximum_frequency, N_freq)
    periods = 1 / frequency

    # Compute the periodogram
    power = ls.power(frequency)

    # Find the index of the peak power
    peak_index = np.argmax(power)

    # Calculate the peak period
    peak_period = periods[peak_index]
    peak_frequency = frequency[peak_index]

    # Find the indices of the top N best periods
    best_indices = np.argsort(power)[-n_best:][::-1]

    # Initialize arrays to store the best periods and values
    best_periods = []
    best_values = []

    # Sort the indices based on the power values in descending order
    sorted_indices = np.argsort(power)[::-1]

    # Iterate over the sorted indices
    for index in sorted_indices:
        current_period = 1 / frequency[index]
        current_value = power[index]

        # Check if the current period is sufficiently different from the
        # existing best periods
        if not any(
            np.abs(current_period - period) / period < periodepsilon
            for period in best_periods
        ):
            best_periods.append(current_period)
            best_values.append(current_value)

        # Break the loop if we have found the desired number of best periods
        if len(best_periods) == n_best:
            break

    # Estimate the noise in the light curve using p2p_rms
    p2p_noise = p2p_rms(fluxes)

    # Evaluate the best-fit sinusoid using the peak period
    best_fit_sinusoid = ls.model(times, peak_frequency)

    # Calculate the model parameters using the peak frequency
    theta = ls.model_parameters(peak_frequency)
    amplitude = theta[1]
    a_90_10_model = (
        np.nanpercentile(best_fit_sinusoid, 90) -
        np.nanpercentile(best_fit_sinusoid, 10)
    )

    # Calculate the chi^2 and reduced chi^2
    residuals = fluxes - best_fit_sinusoid
    chi2 = np.sum((residuals / p2p_noise)**2)
    n_points = len(times)
    n_dof = n_points - 2  # Number of dof (assuming a sinusoid with 2 parameters)
    reduced_chi2 = chi2 / n_dof

    # Create a dictionary to store the results
    # Store the frequency grid and power in a dictionary
    lsp = {
        'times':times, 'fluxs': fluxes,
        'periods': periods, 'power': power, 'bestlspval': np.nanmax(power),
        'nbestperiods': best_periods, 'nbestlspvals': best_values
    }

    baseline = np.nanmax(times) - np.nanmin(times)

    snr_metric = (a_90_10_model / p2p_noise) * np.sqrt(
        baseline / peak_period
    )

    results = {
        'lsp': lsp,
        'period': peak_period,
        'lspval': np.nanmax(power),
        'amplitude': amplitude,
        'a_90_10_model': a_90_10_model,
        't0': np.nanmin(times),
        'outdir': outdir,
        'periodogram_method': periodogram_method,
        'p2p_noise': p2p_noise,
        'best_fit_sinusoid': best_fit_sinusoid,
        'chi2': chi2,
        'reduced_chi2': reduced_chi2,
        'n_points': n_points,
        'n_dof': n_dof,
        'baseline': baseline,
        'snr_metric': snr_metric
    }

    if cachedict is not None:
        for k,v in cachedict.items():
            results[k] = v

    if cache_periodogram_pkls:
        with open(pklpath, 'wb') as f:
            pickle.dump(results, f)
            LOGINFO(f'Made {pklpath}')

    return results
