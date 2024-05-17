import numpy as np, pandas as pd
from os.path import join
from astropy.table import Table
from astropy.io import fits
from copy import deepcopy

localdir = '/Users/luke/local/trp_RESULTS/'
csvpath = join(
    localdir,
    "0to500pc_quicklook_results_trimmed_periodcleaned_supp_nodup_crowd_cut_250pc_cut_20ascrowding.csv"
)
df = pd.read_csv(csvpath)

df['Gaia_DR3_Source_ID'] = 0
df = df.rename({'dr2_source_id': 'Gaia_DR2_Source_ID'}, axis='columns')
df['LegacySurvey_DR8_ID'] = 0
df['PanSTARRS_DR2_ID'] = 0
df['TwoMASS_ID'] = 'NA'
df['delta_ra'] = 0
df['delta_dec'] = 0
df['inertial'] = 0
df['priority'] = 6085
df['cadence'] = 'bright_1x1'
df['instrument'] = 'BOSS'
df['mapper'] = 'MWM'
df['program'] = 'open_fiber'
df['category'] = 'science'
df['cartonname'] = 'manual_mwm_tessperiodic'
df['can_offset'] = 1

selcols = [
    'Gaia_DR3_Source_ID',
    'Gaia_DR2_Source_ID',
    'LegacySurvey_DR8_ID',
    'PanSTARRS_DR2_ID',
    'TwoMASS_ID',
    'ra',
    'dec',
    'delta_ra',
    'delta_dec',
    'inertial',
    'priority',
    'cadence',
    'instrument',
    'mapper',
    'program',
    'category',
    'cartonname',
    'can_offset'
]

sdf1 = df[selcols]

sdf2 = deepcopy(sdf1)
sdf2['instrument'] = 'APOGEE'

df = pd.concat((sdf1, sdf2))

assert df['Gaia_DR2_Source_ID'].isna().sum() == 0

print(len(sdf1))
print(len(df))

# Convert the DataFrame to an Astropy Table
table = Table.from_pandas(df)

# Write the table to a FITS file
fits.writeto('manual_mwm_tessperiodic.fits', table.as_array(), overwrite=True)
