import numpy as np, pandas as pd
from copy import deepcopy
from os.path import join
from cdips.utils.gaiaqueries import given_votable_get_df

# result of job 
#
#   select *
#   from user_lbouma.qlps1s55xgdr2plxgt2 as u, gaiadr3.dr2_neighbourhood as g
#   where u.dr2_source_id = g.dr2_source_id
#
# applied on QLP_s1s55_X_GDR2_parallax_gt_2.csv
#
qlpdir = '/Users/luke/local/QLP'
votablepath = join(qlpdir, 'QLP_dr2_to_dr3_ids-result.vot.gz')

dr2_x_dr3_df = given_votable_get_df(votablepath, assert_equal=None)

df = deepcopy(dr2_x_dr3_df)
del dr2_x_dr3_df

df['abs_magnitude_difference'] = np.abs(df['magnitude_difference'])

get_dr3_xm = lambda _df: (
        _df.sort_values(by='abs_magnitude_difference').
        drop_duplicates(subset='dr2_source_id', keep='first')
)

s_dr3 = get_dr3_xm(df)

outpath = join(
    qlpdir, 'QLP_s1s55_X_GDR2_GDR3_dr2_parallax_gt_2_neighbourhood_result.csv'
)

print(len(s_dr3))

s_dr3.to_csv(outpath, index=False)
print(f"Wrote {outpath}")

import IPython; IPython.embed()
