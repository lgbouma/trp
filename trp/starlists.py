"""
Contents:
    get_ticids
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
import os, pickle, tarfile
from os.path import join
from glob import glob
import numpy as np, pandas as pd

from trp.paths import DATADIR
from trp.tar_utils import extract_tarball

#############

def get_ticids(sample_id, lcpipeline):

    if sample_id == 'debug':
        ticids = [
            "402980664"
        ]

        N_stars_to_search = len(ticids)
        N_lcs_to_search = -1

    elif "iterlist_debug" in sample_id:
        if sample_id == 'iterlist_debug0':
            ticids = [
                "368129164",
                "311092148"
            ]
        elif sample_id == 'iterlist_debug1':
            ticids = [
                "50745567",
                "178155030"
            ]
        else:
            raise NotImplementedError

        N_stars_to_search = len(ticids)
        N_lcs_to_search = -1

    # CONSTRUCTION OF THE JOBS
    elif (
        'pc' in sample_id and 'to' in sample_id
        and lcpipeline=='qlp'
        and 'job' not in sample_id
    ):

        lower = int(sample_id.split("to")[0])
        upper = int(sample_id.split("to")[1].split("pc")[0])

        if lower == 0:
            lower = 1

        assert upper <= 500 # parsecs

        # As described:
        # Generated on wh1 at /ar1/TESS/QLP ;  I downloaded all S1-S55 QLP
        # light curves.  Got the TICID's.  Crossmatched against Gaia DR2 to get
        # the parallaxes.  The TICID's and parallaxes for plx>2
        # (distance<500pc) are inside this file.
        df = pd.read_csv(
            join(DATADIR, "QLP_s1s55_X_GDR2_parallax_gt_2_ticid_plx.tar.gz")
        )

        sel = (
            (df["parallax"] <= 1e3*(1/lower))
            &
            (df["parallax"] > 1e3*(1/upper))
        )

        sdf = df[sel]
        ticids = np.unique(list(sdf["ticid"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = len(sdf) # this is incorrect; in 2024, it's ~3x N_stars_to_search

    elif (
        'pc' in sample_id and 'to' in sample_id
        and lcpipeline=='qlp'
        and 'job' in sample_id
    ):

        # Run on OSG!
        lower = int(sample_id.split("to")[0])
        upper = int(sample_id.split("to")[1].split("pc")[0])

        if lower == 0:
            lower = 1

        assert upper <= 500 # parsecs

        # this file is transferred via the HTCondor submit script
        csvpath = f"./{sample_id}.csv"
        LOGINFO(f"Attempting to get TICIDs from {csvpath}...")
        df = pd.read_csv(csvpath)

        # light curves are passed as tarball via HTCondor..
        tarballpath = f"./{sample_id}_lightcurves.tar.gz"
        extractpath = "./"
        extract_tarball(tarballpath, extractpath)

        ticids = np.unique(list(df["ticid"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = 3*N_stars_to_search # this is very approximate


    # e.g., 30to50pc_mkdwarf, 50to60pc_mkdwarf, etc.
    elif 'pc_mkdwarf' in sample_id and 'to' in sample_id and lcpipeline=='spoc2min':
        raise NotImplementedError

    else:
        raise NotImplementedError

    LOGINFO(42*'-')
    LOGINFO(f"{sample_id}")
    LOGINFO(f"N_stars_to_search = {N_stars_to_search}...")
    LOGINFO(f"N_lcs_to_search = {N_lcs_to_search}...")

    return ticids




