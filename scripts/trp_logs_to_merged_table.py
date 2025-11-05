import os
from trp.pipeline_utils import logpaths_to_csv
import pandas as pd, numpy as np
from glob import glob
from os.path import join
from cdips.utils.gaiaqueries import given_source_ids_get_gaia_data

##########################################
# USAGE: update this
outcsvpath = './20251105_scocenql_v0.csv'
##########################################

logdir = '/ar0/local/trp_cache/rotperiod_finding/unpopular_20251104_scocenql_ticid_list_csv'
logpaths = glob(join(logdir, 'TIC*status.log'))

df = logpaths_to_csv(logpaths, outcsvpath)

dr2_x_tic8_ftrpath = '/home/luke/local/TARS/tic8_plxGT2_TmagLT17_lukebouma.ftr'
xmdf = pd.read_feather(dr2_x_tic8_ftrpath)

mdf = df.merge(xmdf, how='left', left_on='ticid', right_on='ID')
mdf = mdf.rename({'GAIA':'gaia_dr2_source_id'}, axis='columns')

runid = f"dr2_{os.path.basename(outcsvpath).replace('.csv','')}"
dr2_source_ids = np.array([np.int64(mdf.gaia_dr2_source_id)])[0]
gdf = given_source_ids_get_gaia_data(
    dr2_source_ids, runid, n_max=50000, overwrite=True,
    enforce_all_sourceids_viable=True, savstr='', which_columns='*',
    table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=False
)
gdf['source_id'] = gdf['source_id'].astype(str)
gdf = gdf.drop_duplicates(subset='source_id')

mgdf = mdf.merge(gdf, how='left', left_on='gaia_dr2_source_id',
                 right_on='source_id')

selcols = (
    'lcpath,ticid,starid,sector,cadence_sec,periodogram_method,period,amplitude,a_90_10_model,reduced_chi2,bestlspval,nbestperiods,nbestlspvals,p2p_rms,snr_metric,exitcode,gaia_dr2_source_id,Tmag,ra,dec,parallax,pmra,pmdec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,g_rp,radial_velocity,l,b,teff_val'.split(',')
)

outcsvpath = './20251105_scocenql_v0_gdr2_supp.csv'
mgdf[selcols].to_csv(outcsvpath, index=False)
print(f'made {outcsvpath}')

import IPython; IPython.embed()
