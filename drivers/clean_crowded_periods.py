import numpy as np, pandas as pd
from numpy import array as nparr
import matplotlib.pyplot as plt
import os
from os.path import join
from trp.paths import RESULTSDIR

localdir = '/Users/luke/local/trp_RESULTS/'

radepath = join(
    localdir,
    "radecsourceid_0to500pc_quicklook_results_trimmed_periodcleaned_supp_nodup.csv"
)

suppcsvpath = join(
    localdir,
    "0to500pc_quicklook_results_trimmed_periodcleaned_supp_nodup.csv"
)
df = pd.read_csv(suppcsvpath)

if not os.path.exists(radepath):

    selcols = 'ra,dec,dr2_source_id,phot_g_mean_mag'.split(",")
    sdf = df[selcols]
    sdf = sdf.rename({
        'dr2_source_id':'source_id'
    }, axis='columns')

    sdf.to_csv(radepath, index=False)
    print(f'Wrote {radepath}')
    print('Go to the Gaia archive and run the crossmatch.')

# run on gaia archive:
# (a port of given_source_ids_get_neighbor_counts)
#     select top 999999
#     u.source_id, g.source_id, g.ra, g.dec, g.phot_g_mean_mag,
#     3600*DISTANCE(POINT('ICRS', u.ra, u.dec), POINT('ICRS', g.ra, g.dec)) as dist_arcsec,
#     g.phot_g_mean_mag - u.phot_g_mean_mag as d_gmag
#     FROM
#     user_lbouma.trp500pc20240515 as u, gaiadr2.gaia_source as g
#     WHERE
#     1 = CONTAINS(POINT('ICRS', u.ra, u.dec), CIRCLE('ICRS', g.ra, g.dec, 0.002777777777777778))
#     AND
#     g.phot_g_mean_mag - u.phot_g_mean_mag < 2.5
#     AND
#     u.source_id != g.source_id
#     ORDER BY
#     u.source_id, dist_arcsec ASC

source_ids = nparr(df['dr2_source_id']).astype(np.int64)

vot10path = join(localdir, "trp_500pc_nbhrcount_10arcsec-result.vot.gz")
vot20path = join(localdir, "trp_500pc_nbhrcount_20arcsec-result.vot.gz")

from cdips.utils.gaiaqueries import _big_nbhr_count_cleaner

count_df_10, ndf_10 = _big_nbhr_count_cleaner(vot10path, source_ids)
count_df_20, ndf_20 = _big_nbhr_count_cleaner(vot20path, source_ids)

assert len(df) == len(count_df_10)
assert len(df) == len(count_df_20)

df['nbhr_count_10'] = np.array(count_df_10.nbhr_count)
df['flag_dr2_crowding_10'] = (df['nbhr_count_10'] >= 1).astype(int)

df['nbhr_count_20'] = np.array(count_df_20.nbhr_count)
df['flag_dr2_crowding_20'] = (df['nbhr_count_20'] >= 1).astype(int)

suppcsvpath = join(
    localdir,
    "0to500pc_quicklook_results_trimmed_periodcleaned_supp_nodup_crowd.csv"
)
df.to_csv(suppcsvpath, index=False)
print(f"Wrote {suppcsvpath}")

# NOTE if you wanted to remove binaries... could add RUWE....
sel = (
    (df.flag_dr2_crowding_20 == 0)
    &
    ( np.sqrt(df.x**2 + df.y**2 + df.z**2) < 250 )
    &
    ( ~( (df.log10_ampl > -0.25) & (df.period > 7))  )
)
sdf = df[sel]
suppcsvpath = join(
    localdir,
    "0to500pc_quicklook_results_trimmed_periodcleaned_supp_nodup_crowd_cut_250pc_cut_20ascrowding.csv"
)
sdf.to_csv(suppcsvpath, index=False)
print(f"Wrote {suppcsvpath}")

print(sdf.ticid.unique().size)

import IPython; IPython.embed()
