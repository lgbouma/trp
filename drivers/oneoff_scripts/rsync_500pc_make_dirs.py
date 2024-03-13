import os, subprocess
import pandas as pd, numpy as np

def make_dirs():
    csvpath = 's0001_to_s0055_QLP_ticid_path_merged_parallax_gt_2mas.csv'
    df = pd.read_csv(csvpath)

    ticids, paths = np.array(df.ticid), np.array(df.path)
    Ntot = len(ticids)

    for ix, (ticid, relpath) in enumerate(zip(ticids, paths)):

        print(f"{ix}/{Ntot}...")

        dirname = os.path.dirname(relpath)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        # NOTE: initially, the plan was to rsync direct from here.
        # however the relevant OSG singularity container does not have a functioning version of rsync...
        # so the actual rsync is done from bash.

if __name__ == "__main__":
    make_dirs()
