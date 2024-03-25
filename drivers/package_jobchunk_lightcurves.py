"""
Usage: `python package_jobchunk_lightcurves.py`
On: OSG

Iterate over all chunks in a job to create HTCondor-transferrable light curve
tarballs.  (After having made chunklists).
"""

from os.path import join
import os, tarfile
import pandas as pd, numpy as np
from glob import glob
from trp.tar_utils import create_tarball

def package_jobchunk_lightcurves(
    df, volumestr = None,
):

    print(20*'-')
    print(volumestr)

    ENVDATADIR = '/ospool/ap21/data/ekul/chunklists'
    jobdir = join(ENVDATADIR, volumestr)
    assert os.path.exists(jobdir)

    jobpaths = np.sort(glob(join(jobdir, f"{volumestr}_job*_ticid.csv")))
    assert len(jobpaths) >= 1

    for jobpath in jobpaths:

        sample_id = os.path.basename(jobpath).replace(".csv","")
        tarball_path = join(jobdir, f"{sample_id}_lightcurves.tar.gz")
        if os.path.exists(tarball_path):
            print(f'Found {tarball_path}, continue.')
            continue

        _df = pd.read_csv(jobpath)
        mdf = _df.merge(df, on='ticid', how='left')
        mdf = mdf.drop_duplicates(subset=['ticid','path'])

        qlpdir = '/ospool/ap21/data/ekul/QLP/ar1/TESS/QLP'
        fullpaths = [ join(qlpdir, p) for p in mdf.path ]
        fullpathexists = [ os.path.exists(join(qlpdir, p)) for p in mdf.path ]

        mdf['fullpath'] = fullpaths
        mdf['fullpath_exists'] = fullpathexists
        N_exists = np.sum(fullpathexists)

        print(f"{sample_id}: {len(_df)} stars, "
              f"{len(mdf)} LCs expected, "
              f"{N_exists} exist.")

        create_tarball(fullpaths, tarball_path)


if __name__ == "__main__":

    # dataframe of ticid's and (relative) light curve paths
    df = pd.read_csv(
        "/ospool/ap21/data/ekul/QLP/"
        "s0001_to_s0055_QLP_ticid_path_merged_parallax_gt_2mas.csv"
    )

    #upper_grid = np.arange(5, 505, 5)
    #upper_grid = np.arange(5, 25, 5)
    lower_grid = np.arange(0, 500, 5)
    upper_grid = np.arange(5, 105, 5)

    lower_grid = np.arange(100, 500, 5)
    upper_grid = np.arange(105, 505, 5)


    for lower, upper in zip(lower_grid, upper_grid):
        volumestr = f'{lower}to{upper}pc'
        package_jobchunk_lightcurves(df, volumestr=volumestr)
