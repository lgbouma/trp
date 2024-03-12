"""
Usage: `python make_chunklist.py`
On: mico

Generate input lists needed to send thousands of jobs (each with a manageable
number of stars) to condor.

Each job contains a subset of a volume slice; at most `chunksize` stars.

Volume slices are by successive 5 parsec distance cuts, to avoid overwhelming
the condor queue.

(These lists chunk ticid's from QLP)
"""

from trp.starlists import get_ticids
from trp.paths import DATADIR
from os.path import join
import os
import pandas as pd, numpy as np
from glob import glob

def chunk_list(lst, chunksize):
    return [lst[i:i+chunksize] for i in range(0, len(lst), chunksize)]

def make_job_chunks(
    volumestr = None,
    chunksize = 50,
    max_n_jobs = 9999
):

    assert isinstance(volumestr, str)

    outdir = "/Users/luke/local/trp/chunklists"
    if not os.path.exists(outdir): os.mkdir(outdir)
    outdir = join(outdir, volumestr)
    if os.path.exists(outdir):
        print(f"Found {outdir}")
        return outdir

    if not os.path.exists(outdir): os.mkdir(outdir)

    ticids = get_ticids(volumestr, 'qlp')
    print(f'{volumestr}: Got {len(ticids)} TICIDs')

    # can have at most ~50k stars per volume slice given current OSG job number
    # limits...
    assert len(ticids) <= max_n_jobs * chunksize

    # list of "chunks"; each chunk represents one condor job (and one
    # "sampleid" in the trp_pipeline verbiage.  max ~50 stars per chunk to fit
    # within 10 hour wallclock time max per job.
    chunked_ticids = chunk_list(ticids, chunksize)

    assert len(chunked_ticids) <= max_n_jobs

    for job_ix, chunk in enumerate(chunked_ticids):

        outcsv = join(
            outdir, f"{volumestr}_job{str(job_ix).zfill(7)}_ticid.csv"
        )
        df = pd.DataFrame({
            'ticid': chunk
        })
        df.to_csv(outcsv, index=False)
        print(f"Wrote {outcsv}")

    return outdir

def make_job_lists(outdir, volumestr):

    basedir = os.path.dirname(outdir)
    csvpaths = np.sort(glob(join(outdir, "*.csv")))
    namelist = [os.path.basename(c.rstrip(".csv")) for c in csvpaths]

    outpath = join(basedir, volumestr+".sampleids")

    with open(outpath, 'w') as f:
        f.writelines(f"{name}\n" for name in namelist)

    print(f'Wrote {outpath}')


if __name__ == "__main__":

    lower_grid = np.arange(0, 500, 5)
    upper_grid = np.arange(5, 505, 5)

    for lower, upper in zip(lower_grid, upper_grid):

        volumestr = f'{lower}to{upper}pc'

        outdir = make_job_chunks(volumestr=volumestr, chunksize=50)

        make_job_lists(outdir, volumestr)

