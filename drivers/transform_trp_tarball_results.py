"""
Given a directory with a bunch of joboutput_130to135pc_job0000185_ticid.tar.gz
style files generated by the trp pipeline... unpackage them, read the output to
a CSV cache file, and clean.
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

from glob import glob
from os.path import join
import os
import pandas as pd, numpy as np
from collections import Counter

from trp.tar_utils import extract_tarball
import trp.pipeline_utils as pu

def build_pipeline_dataframe(tarballdir, do_tarball_extraction=1):
    """
    Given a directory with a bunch of
    joboutput_130to135pc_job0000185_ticid.tar.gz style files generated by the
    trp pipeline... unpackage them, read the output to a CSV cache file.
    """

    tarballpaths = np.sort(glob(join(tarballdir, "*.tar.gz")))

    if len(tarballpaths) > 0 and do_tarball_extraction:
        for tarballpath in tarballpaths:
            extract_tarball(tarballpath, tarballdir)

    logpaths = np.sort(glob(join(tarballdir, "srv", "joboutput*", "*log")))

    outdir = os.path.dirname(tarballdir)
    outcsvpath = join(outdir, "0to500pc_quicklook_results.csv")
    trimoutcsvpath = join(outdir, "0to500pc_quicklook_results_trimmed.csv")

    if not os.path.exists(trimoutcsvpath):
        colnames = (
            'lcpath,ticid,starid,sector,cadence_sec,periodogram_method,'
            'period,a_90_10_model,bestlspval,reduced_chi2,t0,nbestperiods,'
            'nbestlspvals,exitcode'.split(",")
        )

        df = pu.logpaths_to_csv(logpaths, outcsvpath, colnames=colnames)
        df['log10_ampl'] = np.log10(df['a_90_10_model'])
        df.to_csv(outcsvpath, index=False)
        LOGINFO(f"Wrote {outcsvpath}")

        selcols = (
            'ticid,period,bestlspval,a_90_10_model,reduced_chi2,'
            'sector,log10_ampl'.split(",")
        )
        sdf = df[selcols]
        sdf.to_csv(trimoutcsvpath, index=False)
        LOGINFO(f"Wrote {trimoutcsvpath}")
    else:
        sdf = pd.read_csv(trimoutcsvpath)

    # Count the number of entries with different exitcode values
    exitcode_counts = Counter(df['exitcode'])
    LOGINFO("Exitcode counts:")
    for exitcode, count in exitcode_counts.items():
        LOGINFO(f"Exitcode {exitcode}: {count} entries")




if __name__ == "__main__":

    ##########################################
    # USER OPTIONS HERE #
    tarballdir = '/Users/luke/local/trp_RESULTS/RESULTS'
    do_tarball_extraction = 0
    trimoutcsvpath = join(
        os.path.dirname(tarballdir),
        "0to500pc_quicklook_results_trimmed.csv"
    )
    ##########################################

    if not os.path.exists(trimoutcsvpath):
        build_pipeline_dataframe(
            tarballdir, do_tarball_extraction=do_tarball_extraction
        )

    df = pd.read_csv(trimoutcsvpath)
    if 'log10_ampl' not in df.columns:
        df['log10_ampl'] = np.log10(df['a_90_10_model'])

    gcsvpath = '/Users/luke/local/QLP/QLP_s1s55_X_GDR2_parallax_gt_2.csv'
    gdf = pd.read_csv(gcsvpath)
    selcols = 'ra,dec,parallax,pmra,pmdec,M_G,bp_rp,ticid'.split(",")
    gdf = gdf[selcols]

    mdf = df.merge(gdf, on='ticid', how='inner')

    from sunnyhills.physicalpositions import calculate_XYZ_given_RADECPLX
    x,y,z = calculate_XYZ_given_RADECPLX(mdf.ra, mdf.dec, mdf.parallax)
    mdf['x'] = x
    mdf['y'] = y
    mdf['z'] = z

    suppcsvpath = join(os.path.dirname(tarballdir),
                       "0to500pc_quicklook_results_trimmed_supp.csv")
    mdf.to_csv(suppcsvpath, index=False)
    print(f"Wrote {suppcsvpath}")

    mdf = mdf.sort_values(by='bestlspval', ascending=True)
    mdf = mdf.drop_duplicates(subset='ticid', keep='first')

    suppcsvpath = join(os.path.dirname(tarballdir),
                       "0to500pc_quicklook_results_trimmed_supp_nodup.csv")
    mdf.to_csv(suppcsvpath, index=False)
    print(f"Wrote {suppcsvpath}")
