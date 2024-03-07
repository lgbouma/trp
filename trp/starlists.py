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
import os, pickle
from os.path import join
from glob import glob
import numpy as np, pandas as pd

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

    elif 'pc' in sample_id and 'to' in sample_id and lcpipeline=='qlp':

        lower = int(sample_id.split("to")[0])
        upper = int(sample_id.split("to")[1].split("pc")[0])

        assert upper <= 500 # parsecs

        df = pd.read_csv(join(QLPDIR, "QLP_s1s55_X_GDR2_parallax_gt_2.csv"))

        sel = (
            (df["parallax"] <= 1e3*(1/lower))
            &
            (df["parallax"] > 1e3*(1/upper))
        )

        sdf = df[sel]
        ticids = np.unique(list(sdf["ticid"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = len(sdf) # this is incorrect; in 2023, it's ~3x N_stars_to_search

    # e.g., 30to50pc_mkdwarf, 50to60pc_mkdwarf, etc.
    elif 'pc_mkdwarf' in sample_id and 'to' in sample_id and lcpipeline=='spoc2min':

        lower = int(sample_id.split("to")[0])
        upper = int(sample_id.split("to")[1].split("pc")[0])

        df = pd.read_csv(join(SPOCDIR, "gaia_X_spoc2min_merge.csv"))

        sel = (
            (df["M_G"] > 4)
            &
            (df["bp_rp"] > 1.5)
            &
            (df["TESSMAG"] < 16)
            &
            (df["parallax"] <= 1e3*(1/lower))
            &
            (df["parallax"] > 1e3*(1/upper))
        )

        sdf = df[sel]
        ticids = np.unique(list(sdf["TICID"].astype(str)))

        N_stars_to_search = len(ticids)
        N_lcs_to_search = len(sdf)

    else:
        raise NotImplementedError

    LOGINFO(42*'-')
    LOGINFO(f"{sample_id}")
    LOGINFO(f"N_stars_to_search = {N_stars_to_search}...")
    LOGINFO(f"N_lcs_to_search = {N_lcs_to_search}...")

    return ticids




