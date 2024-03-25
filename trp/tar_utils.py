"""
Utilities for zipping and unzipping tarballs.

Contents:
    create_tarball
    extract_tarball
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
from os.path import join
import os, tarfile
import pandas as pd, numpy as np
from glob import glob

def create_tarball(fullpaths, tarball_path):
    """
    Given a list of paths, write them into a gzipped tar achive.
    """
    with tarfile.open(tarball_path, "w:gz") as tar:
        for file_path in fullpaths:
            if os.path.isfile(file_path):
                tar.add(file_path, arcname=os.path.basename(file_path))
            else:
                LOGINFO(f"Warning: {file_path} is not a valid file.")
    LOGINFO(f"...Made {tarball_path}")


def extract_tarball(tarball_name, extract_path):
    """
    Unzip a gzipped tar archive.
    """
    with tarfile.open(tarball_name, "r:gz") as tar:
        tar.extractall(path=extract_path)
        LOGINFO(f"Extracted {tarball_name} to {extract_path}")
