"""
Contents:
    run_trp
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

#TODO
#from complexrotators.lcprocessing import (
#    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
#)
#from complexrotators.plotting import (
#    plot_phased_light_curve, plot_dipcountercheck, plot_cpvvetter
#)


#############

def run_trp():

    #################
    # begin options #
    #################
    forcepdf = 1 # if yes, perhaps also have "LOCALDEBUG" set true..

    lcpipeline = 'qlp' # "qlp" or "spoc2min"

    sample_ids = [
        'debug'
        # ### example samples:
        #'1to20pc'
        # '20to40pc'
    ]

    ###############
    # end options #
    ###############

    for sample_id in sample_ids:

        ticids = get_ticids(sample_id, lcpipeline)

        if len(ticids) > 100 and forcepdf:
            raise NotImplementedError

        for ticid in ticids:
            LOGINFO(42*'-')
            LOGINFO(f"Beginning {ticid}...")
            find_rotperiod(
                ticid, sample_id, forcepdf=forcepdf, lcpipeline=lcpipeline
            )

    LOGINFO("Finished 🎉🎉🎉")


def find_rotperiod(ticid, sample_id, forcepdf=0, lcpipeline='qlp'):
    """
    Args:

    ticid: e.g. "289840928"

    sample_id: e.g., "30to50pc_mkdwarf" (used for cacheing)

    forcepdf: if true, will require the pdf plot to be made, even if the usual
        exit code criteria were not met.

    lcpipeline: "qlp" or "spoc2min"

    exit code definitions:
        exitcode 1: periodogram_condition was met.  e.g.,
            periodogram_condition = (period < 10) & (lspval > 0.1)

        exitcode 2: periodogram_condition was not met.

        exitcode 3: (todo)

        exitcode 4: (todo)
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

    MINIMUM_EXITCODE = 2
    if forcepdf:
        MINIMUM_EXITCODE = 1
    #1 if any kind of exit means do not rerun
    #2 if only a periodogram or not enoigh dip exit means dont rerun
    if minexitcode >= MINIMUM_EXITCODE:
        LOGINFO(f"TIC{ticid}: found log for {ticid} with exitcode {minexitcode}. skip.")
        return 0


    #
    # get data
    #
    lcpaths = get_lcpaths_fromlightkurve_given_ticid(ticid, lcpipeline)

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

        #import IPython; IPython.embed()
        # #FIXME TODO: implement!
        # # get the relevant light curve data
        # (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
        #  sector, starid) = prepare_rot_light_curve(
        #      lcpath, cachedir, lcpipeline=lcpipeline
        #  )


