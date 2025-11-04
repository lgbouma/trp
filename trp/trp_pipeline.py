"""
Contents:
    run_trp
    find_rotperiod
    prepare_rot_light_curve

Minor helpers:
    parse_sample_id
"""
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
import os, pickle, re
from os.path import join
from glob import glob
import numpy as np, pandas as pd

from trp.paths import DATADIR, RESULTSDIR, CACHEDIR
from trp import pipeline_utils as pu

from trp.starlists import get_ticids
from trp.getters import get_lcpaths_fromlightkurve_given_ticid
from trp.lcprocessing import (
    rotation_periodsearch, calculate_lsp, time_bin_lightcurve
)

from complexrotators.lcprocessing import prepare_cpv_light_curve
from astropy.io import fits

AESTHETIC_IMPORT_WORKS = 0
try:
    from aesthetic.plot import savefig, format_ax, set_style
    AESTHETIC_IMPORT_WORKS = 1
except:
    LOGINFO('Failed to import aesthetic; skipping.')

#############

def run_trp(sample_id, mask_known_transits=False, write_vetplot=False,
            forcerun=False, lcpipeline='qlp', periodogram_method='astropyls'):
    """
    lcpipeline:
        'qlp', 'spoc2min', or 'unpopular'.  If 'unpopular', need to have
        created your light curves separately.  'qlp' and 'spoc2min' has logic
        to look on harddrive.

    periodogram_method:
        'astropyls' is the only fully-implemented option

    sample_id:
        This unique identifying string pairs to a list including at least one
        ticid, via `trp.starlists.get_ticids`.

    Implemented sample_ids:
        * 'iterlist_debug0',
        * 'iterlist_debug1',
        * 'debug'
        * "NtoMpc" (e.g., '1to10pc', '10to20pc', etc)
            (for scale, 90to100pc is ~35k stars, ~100k lcs)
    """

    #################
    # begin options #
    #################
    cache_periodogram_pkls = 1

    ###############
    # end options #
    ###############

    sample_id = parse_sample_id(sample_id)

    ticids = get_ticids(sample_id, lcpipeline)

    if len(ticids) > 100 and forcerun:
        raise NotImplementedError

    for ticid in ticids:
        LOGINFO(42*'-')
        LOGINFO(f"Beginning {ticid}...")
        find_rotperiod(
            ticid, sample_id, forcerun=forcerun, lcpipeline=lcpipeline,
            periodogram_method=periodogram_method,
            write_vetplot=write_vetplot,
            cache_periodogram_pkls=cache_periodogram_pkls,
            mask_known_transits=mask_known_transits
        )

    LOGINFO("Finished 🎉🎉🎉")


def parse_sample_id(sample_id):
    """
    If a CSV file of TICIDs is given, ensure formatting is correct.
    """

    if sample_id.endswith(".csv"):
        sample_id = sample_id.replace(".", "_")

        from trp.paths import TARGETDIR
        csvpath = join(TARGETDIR, sample_id.replace("_csv", ".csv"))

        assert os.path.exists(csvpath), f"Could not find {csvpath}"

        msg = f"'ticid' not found in the first line of {csvpath}"
        with open(csvpath, 'r') as file:
            assert 'ticid' in file.readline().strip().split(','), msg

    return sample_id


def find_rotperiod(ticid, sample_id, forcerun=0, lcpipeline='qlp',
                   periodogram_method='astropyls', write_vetplot=0,
                   cache_periodogram_pkls=1, mask_known_transits=0):
    """
    This pipeline takes a light curve (SPOC 2-minute, QLP, or unpopular),
    removes non-zero quality flags, and median-normalizes.  It then bins to a
    30 minute cadence, and runs a periodogram (by default, lomb scargle) on the
    resulting points.  Results are cached in a mix of text files using
    `configparser` and pickle files.

    Args:

        ticid (str): e.g. "289840928"

        sample_id (str): e.g., "30to50pc_mkdwarf" (used for cacheing)

        forcerun (bool): if true, will require the pipeline to be run, even if
            the usual exit code criteria were not met.

        lcpipeline (str): "qlp" or "spoc2min"

        periodogram_method (str): "ls", "astropyls", or "pdm"

        write_vetplot (bool): whether to generate a PNG that can be used to
            vet the validity for a proposed rotation period.

        cache_periodogram_pkls (bool): save pickle files with detailed
        periodogram info.

    exit code definitions:

        exitcode 0: a periodgram was successfully calculated.

        exitcode 1: failed; non-finite light curve.

        exitcode 2: failed; too few points for a period.

        exitcode 3: failed to allocate memory
    """

    if write_vetplot:
        msg = 'You need to install https://github.com/lgbouma/aesthetic'
        assert AESTHETIC_IMPORT_WORKS, msg

    #
    # set up / parse log files
    #

    assert lcpipeline in ["qlp", "spoc2min", "unpopular"]

    cachedir = join(CACHEDIR, "rotperiod_finding")
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    cachename = f"{lcpipeline}_{sample_id}"
    cachedir = join(cachedir, cachename)
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    minexitcode = -1
    cand_logpaths = glob(join(cachedir, f"*tess*00{ticid}*runstatus.log"))
    foundexitcodes = []
    if len(cand_logpaths) > 0:
        for cand_logpath in cand_logpaths:
            st = pu.load_status(cand_logpath)
            if 'exitcode' in st:
                exitcode = st['exitcode']['exitcode']
                foundexitcodes.append(int(exitcode))
        if len(foundexitcodes) > 0:
            minexitcode = np.nanmin(foundexitcodes)

    MINIMUM_EXITCODE = 0
    if forcerun:
        MINIMUM_EXITCODE = -1
    if minexitcode >= MINIMUM_EXITCODE and not forcerun:
        LOGINFO(f"TIC{ticid}: found log {cand_logpath} for {ticid} with "
                f"exitcode {minexitcode}. skip.")
        return 0

    #
    # get data
    #
    lcpaths = []
    if lcpipeline == 'qlp':
        # Assuming HTCondor transferred tarball, and extraction happened in
        # get_ticids
        LOGINFO(f"Beginning local fileglob search for TIC {ticid}...")
        lcpaths = glob(f"./hlsp_qlp_tess_ffi_s*-0*{ticid}_tess*llc.fits")
        LOGINFO(f"Got N={len(lcpaths)} for TIC {ticid}...")

    if len(lcpaths) == 0 and lcpipeline == 'qlp':
        # Fall back to lightkurve attempt if transfer failed.
        lcpaths = get_lcpaths_fromlightkurve_given_ticid(
            ticid, lcpipeline, cachedir=cachedir, require_lc=0
        )

    if lcpipeline == 'unpopular':
        LOGINFO(f"Beginning unpopular fileglob search for TIC {ticid}...")
        # a hack, for now
        LCDIR = '/Users/luke/Dropbox/proj/wrapunpopular/results/sco-cen-quicklook'
        lcpaths = glob(join(LCDIR, f"*{ticid}*_cpm_llc.csv"))
        LOGINFO(f"Got N={len(lcpaths)} for TIC {ticid}...")

    if len(lcpaths) == 0:
        LOGINFO(f"TIC{ticid}: Failed to get any light curves; continue.")
        return 0

    #
    # for each light curve (sector / cadence specific), clean, calculate the
    # periodograms, and save output.
    #
    for lcpath in lcpaths:

        assert os.path.exists(lcpath)

        # instantiate the log
        lcpbase = re.sub(r'\.(fits|csv)$', '', os.path.basename(lcpath))
        logpath = join(cachedir, f'{lcpbase}_runstatus.log')
        if not os.path.exists(logpath):
            lcpd = {
                'lcpath': lcpath,
                'ticid': ticid
            }
            pu.save_status(logpath, 'lcpath', lcpd)
            LOGINFO(f"Made {logpath}")

        st = pu.load_status(logpath)
        if 'exitcode' in st:
            exitcode = st['exitcode']['exitcode']
            if minexitcode >= MINIMUM_EXITCODE:
                LOGINFO(f"{lcpbase}: found exitcode {exitcode}. skip.")
                if not forcerun:
                    continue

        if lcpath.endswith(".fits"):
            with fits.open(lcpath) as hdulist:
                hdr = hdulist[0].header
            selcols = "RA_OBJ,DEC_OBJ,TESSMAG,RADIUS,LOGG,MASS,TEFF".split(",")
            starinfod = {k:hdr[k] for k in selcols}
            pu.save_status(logpath, 'starinfo', starinfod)

        (time, flux, qual, x_obs, y_obs, _, _, _, cadence_sec,
         sector, starid) = prepare_rot_light_curve(
             lcpath, cachedir, lcpipeline=lcpipeline
         )

        if mask_known_transits:
            # If the system has known planet(s), mask the transits.
            # TODO: implement
            raise NotImplementedError('Need to finish implementing this')
            from trp.getters import get_ephemeris_given_ticid
            ephem_df = get_ephemeris_given_ticid(ticid)
            #df.columns = ["Epoch", "Period", "Duration"]
            t0 = ephem_df["Epoch"].iloc[0]
            period = ephem_df["Period"].iloc[0]
            dur = ephem_df["Duration"].iloc[0]  / 24

        if y_obs is None:
            LOGWARNING(f"{starid}: Failed; non-finite light curve.")
            exitcode = {'exitcode': 1}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue

        bin_size_days = 30 / (60*24)
        btime, bflux = time_bin_lightcurve(x_obs, y_obs, bin_size_days)

        if np.sum(np.isfinite(y_obs)) <= 20:
            LOGWARNING(f"{starid}: Failed; too few points for a period.")
            exitcode = {'exitcode': 2}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue

        #
        # Calculate the periodogram using the binned (30-minute cadence) light
        # curve.  Cache the original (non-zero quality flags, and
        # median-normalized) light curve for future viz, x_obs and y_obs.
        #
        nworkers = 1 # 1 # if "None" will multithread
        cachedict = {
            'x_obs': x_obs, # as above.
            'y_obs': y_obs
        }
        try:
            if periodogram_method == 'astropyls':
                d = calculate_lsp(
                    btime, bflux, starid, cachedir, cachedict=cachedict
                )
            if periodogram_method == 'ls':
                raise NotImplementedError('will break amplitudes below')
                d = rotation_periodsearch(
                    btime, bflux, starid, cachedir, t0='binmin',
                    periodogram_method=periodogram_method, do_finetune=0,
                    write_pngs=0,
                    nworkers=nworkers, cachedict=cachedict
                )
        except OSError as e:
            LOGWARNING(f"{starid}: Failed to allocate memory ({e}).")
            exitcode = {'exitcode': 3}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue

        psr = {
            'starid': starid,
            'sector': sector,
            'cadence_sec': cadence_sec,
            'periodogram_method': periodogram_method,
            'period': d['period'], # default lomb scargle period
            'amplitude': np.abs(d['amplitude']),
            'a_90_10_model': d['a_90_10_model'],
            'reduced_chi2': d['reduced_chi2'],
            'bestlspval': d['lsp']['bestlspval'],
            't0': d['t0'],
            'nbestperiods': d['lsp']['nbestperiods'],
            'nbestlspvals': d['lsp']['nbestlspvals'],
            'p2p_rms': d['p2p_noise'],
            'snr_metric': d['snr_metric'],
        }
        pu.save_status(logpath, f'rotation_periodsearch_{periodogram_method}_results', psr)
        LOGINFO(f"Updated {logpath} with "
                f"rotation_periodsearch_{periodogram_method}_results")

        LOGINFO(f"{starid}: saved {periodogram_method} results; finished.")
        exitcode = {'exitcode': 0}
        pu.save_status(logpath, 'exitcode', exitcode)

        if write_vetplot:
            outpath = join(cachedir, f'{starid}_rotvetter.pdf')
            if not os.path.exists(outpath):
                from trp.plotting import plot_rotvetter
                plot_rotvetter(
                    outpath, lcpath, starid, periodsearch_result=d,
                    lcpipeline=lcpipeline
                )
            else:
                LOGINFO(f"Found {outpath}")


def prepare_rot_light_curve(lcpath, cachedir, lcpipeline='qlp'):
    """
    Given a light curve (SPOC 2-minute or QLP), remove non-zero quality flags,
    and median-normalize.
    Cache the output.
    """

    return prepare_cpv_light_curve(
        lcpath, cachedir, returncadenceno=0, lcpipeline=lcpipeline,
        runmedianfilter=0, rotmode=1
    )
