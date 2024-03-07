"""
Contents:
    get_lcpaths_fromlightkurve_given_ticid
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
import lightkurve as lk

from trp.paths import CACHEDIR

#############

def get_lcpaths_fromlightkurve_given_ticid(
    ticid, lcpipeline, require_lc=1, cachedir=CACHEDIR
):
    """
    Args:
        ticid: e.g. "289840928"
        lcpipeline: "qlp", "spoc2min", or "cdips'

    Returns:
        list of light curve paths
    """
    # ticid like "289840928"

    assert isinstance(ticid, str)
    if ticid.startswith("TIC"):
        ticid_str = ticid
    else:
        ticid_str = f"TIC {ticid}"

    LOGINFO(f"Beginning search for {ticid_str}...")
    lcset = lk.search_lightcurve(ticid_str)

    if lcpipeline == 'spoc2min':
        sel = (lcset.author=='SPOC') & (lcset.exptime.value == 120)
    elif lcpipeline == 'qlp':
        sel = (lcset.author=='QLP')
    elif lcpipeline == 'cdips':
        sel = (lcset.author=='CDIPS')
    lcc = lcset[sel].download_all(download_dir=cachedir)

    if lcpipeline == 'spoc2min':
        lcpaths = glob(
            join(cachedir, f'tess*{ticid}*-s', f'tess*{ticid}*-s_lc.fits')
        )
    elif lcpipeline in ['qlp', 'cdips']:
        lcpaths = glob(
            join(cachedir, "mast*", "HLSP", f'hlsp_{lcpipeline}*{ticid}*', f'*{ticid}*.fits')
        )

    if require_lc:
        msg = f'{ticid}: did not get LC'
        assert len(lcpaths) > 0, msg

    LOGINFO(f"Returning {len(lcpaths)} light curves for {ticid_str}...")

    return lcpaths
