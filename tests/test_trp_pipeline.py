import os
from glob import glob
from os.path import join

###############
# CLEAR CACHE #
###############
from trp.paths import CACHEDIR

lcpipeline = 'qlp'
sample_id = 'debug_test_0'

cachedir = join(CACHEDIR, "rotperiod_finding")
cachename = f"{lcpipeline}_{sample_id}"
cachedir = join(cachedir, cachename)
cand_logpaths = glob(join(cachedir, f"*tess*00*.log"))

if len(cand_logpaths) > 0:
    for cand_logpath in cand_logpaths:
        os.remove(cand_logpath)

################
# RUN PIPELINE #
################
from trp.trp_pipeline import run_trp

run_trp('debug_test_0')

cand_logpaths = glob(join(cachedir, f"*tess*s0044*00*.log"))
assert len(cand_logpaths) == 1

from trp.pipeline_utils import status_to_dict

d = status_to_dict(cand_logpaths[0])

assert abs(float(d['period']) - 2.95) < 0.01
