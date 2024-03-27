"""
Contents:
    save_status
    load_status
    status_to_dict
    logpaths_to_csv
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
import os
from os.path import join
from glob import glob
import configparser
import numpy as np, pandas as pd

def save_status(status_file, section, state_vars):
    """
    Save pipeline status

    Args:

        status_file (string): name of output file

        section (string): name of section to write

        state_vars (dict): dictionary of all options to populate the specified
        section
    """

    config = configparser.RawConfigParser()

    if os.path.isfile(status_file):
        config.read(status_file)

    if not config.has_section(section):
        config.add_section(section)

    for key, val in state_vars.items():
        config.set(section, key, val)

    with open(status_file, 'w') as f:
        config.write(f)


def load_status(status_file):
    """
    Load pipeline status

    Args:
        status_file (string): name of configparser file

    Returns:
        configparser.RawConfigParser
    """

    config = configparser.RawConfigParser()
    gl = config.read(status_file)

    return config


def status_to_dict(
    status_file,
    colnames = (
        'lcpath,ticid,starid,sector,cadence_sec,periodogram_method,'
        'period,bestlspval,t0,nbestperiods,nbestlspvals,exitcode'.split(",")
    )
):

    st = load_status(status_file)

    outdict = {}
    for colname in colnames:
        outdict[colname] = np.nan

    sections = st.sections()
    for section in sections:
        keys = list(st[section].keys())
        for key in keys:
            val = st[section][key]
            if key in colnames:
                outdict[key] = val

    return outdict


def chunk_list(lst, chunksize):
    return [lst[i:i+chunksize] for i in range(0, len(lst), chunksize)]


def logpaths_to_csv(
    logpaths,
    outcsvpath,
    colnames = (
        'lcpath,ticid,starid,sector,cadence_sec,periodogram_method,'
        'period,bestlspval,t0,nbestperiods,nbestlspvals,exitcode'.split(",")
    ),
    chunksize=1000
):

    cachedir = os.path.dirname(outcsvpath)
    if not os.path.exists(cachedir): os.mkdir(cachedir)
    cachedir = join(cachedir, "temp_chunkcsvs")
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    N = len(logpaths)
    LOGINFO(f"Got N={N} logpaths...")

    chunked_logpaths = chunk_list(logpaths, chunksize)
    N_chunks = len(chunked_logpaths)
    LOGINFO(f"-> N={N_chunks} chunks w/ chunksize={chunksize}...")

    for ix, chunked_logpath_list in enumerate(chunked_logpaths):
        rows = []
        for logpath in chunked_logpath_list:
            row = status_to_dict(logpath, colnames=colnames)
            rows.append(row)
        df = pd.DataFrame(rows)
        outpath = join(cachedir, f"chunk_{str(ix).zfill(8)}.csv")
        df.to_csv(outpath, index=False)
        LOGINFO(f"Chunk {ix}/{N_chunks}: made {outpath}")

    tempcsvpaths = np.sort(glob(join(cachedir, "chunk_*csv")))
    df = pd.concat((pd.read_csv(f) for f in tempcsvpaths))
    return df
