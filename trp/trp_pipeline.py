"""
Contents:
    run_trp
    find_rotperiod
    prepare_rot_light_curve
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
import os, pickle
from os.path import join
from glob import glob
import numpy as np, pandas as pd

from trp.paths import DATADIR, RESULTSDIR, CACHEDIR
from trp import pipeline_utils as pu

from trp.starlists import get_ticids
from trp.getters import get_lcpaths_fromlightkurve_given_ticid
from trp.lcprocessing import rotation_periodsearch

from complexrotators.lcprocessing import prepare_cpv_light_curve
from astrobase.lcmath import time_bin_magseries

#############

def run_trp(sample_id):
    """
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
    forcepdf = 0 # if yes, perhaps also have "LOCALDEBUG" set true..
    write_astrobase_pngs = 0

    lcpipeline = 'qlp' # "qlp" or "spoc2min"
    periodogram_method = 'ls'
    ###############
    # end options #
    ###############

    ticids = get_ticids(sample_id, lcpipeline)

    if len(ticids) > 100 and forcepdf:
        raise NotImplementedError

    for ticid in ticids:
        LOGINFO(42*'-')
        LOGINFO(f"Beginning {ticid}...")
        find_rotperiod(
            ticid, sample_id, forcepdf=forcepdf, lcpipeline=lcpipeline,
            periodogram_method=periodogram_method,
            write_astrobase_pngs=write_astrobase_pngs
        )

    LOGINFO("Finished 🎉🎉🎉")


def find_rotperiod(ticid, sample_id, forcepdf=0, lcpipeline='qlp',
                   periodogram_method='ls', write_astrobase_pngs=0):
    """
    This pipeline takes a light curve (SPOC 2-minute or QLP), remove non-zero
    quality flags, and median-normalizes.  It then bins to a 30 minute cadence,
    and runs a periodogram (by default, lomb scargle) on the resulting points.
    Results are cached in a mix of text files using `configparser` and pickle
    files.

    Args:

        ticid (str): e.g. "289840928"

        sample_id (str): e.g., "30to50pc_mkdwarf" (used for cacheing)

        forcepdf (bool): if true, will require the pdf plot to be made, even if the usual
            exit code criteria were not met.

        lcpipeline (str): "qlp" or "spoc2min"

        periodogram_method (str): "ls" or "pdm"

        write_astrobase_pngs (bool): whther to generate the astrobase checkplot pngs,
        a 3x3 grid showing the light curve, periodogram, and phased versions of
        the light curve.  Good for debugging; not good enough for assessing
        what is really happening.

    exit code definitions:

        exitcode 0: a periodgram was successfully calculated.

        exitcode 1: failed; non-finite light curve.

        exitcode 2: failed; too few points for a period.

        exitcode 3: failed to allocate memory
    """

    #
    # set up / parse log files
    #

    assert lcpipeline in ["qlp", "spoc2min"]

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
    if forcepdf:
        MINIMUM_EXITCODE = -1
    if minexitcode >= MINIMUM_EXITCODE:
        LOGINFO(f"TIC{ticid}: found log for {ticid} with exitcode {minexitcode}. skip.")
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
        lcpbase = os.path.basename(lcpath).replace(".fits", "")
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
                if not forcepdf:
                    continue

        (time, flux, qual, x_obs, y_obs, _, _, _, cadence_sec,
         sector, starid) = prepare_rot_light_curve(
             lcpath, cachedir, lcpipeline=lcpipeline
         )

        if y_obs is None:
            LOGWARNING(f"{starid}: Failed; non-finite light curve.")
            exitcode = {'exitcode': 1}
            pu.save_status(logpath, 'exitcode', exitcode)
            continue

        bd = time_bin_magseries(x_obs, y_obs, minbinelems=1)
        btime, bflux = bd['binnedtimes'], bd['binnedmags']

        if bd['nbins'] <= 20:
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
            d = rotation_periodsearch(
                btime, bflux, starid, cachedir, t0='binmin',
                periodogram_method=periodogram_method, do_finetune=0,
                write_pngs=write_astrobase_pngs,
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
            'bestlspval': d['lsp']['bestlspval'],
            't0': d['t0'],
            'nbestperiods': d['lsp']['nbestperiods'],
            'nbestlspvals': d['lsp']['nbestlspvals']
        }
        pu.save_status(logpath, f'rotation_periodsearch_{periodogram_method}_results', psr)
        LOGINFO(f"Updated {logpath} with "
                f"rotation_periodsearch_{periodogram_method}_results")

        LOGINFO(f"{starid}: saved {periodogram_method} results; finished.")
        exitcode = {'exitcode': 0}
        pu.save_status(logpath, 'exitcode', exitcode)

        ##########################################

        # NOTE: may wish to add option to generate a vetting plot in the
        # future; for now, omit and just calculate periodograms.

        # outpath = join(cachedir, f'{starid}_cpvvetter.pdf')
        # if not os.path.exists(outpath):
        #     plot_cpvvetter(
        #         outpath, lcpath, starid, periodsearch_result=d,
        #         findpeaks_result=r, lcpipeline=lcpipeline
        #     )
        # else:
        #     LOGINFO(f"Found {outpath}")

        ##########################################


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
