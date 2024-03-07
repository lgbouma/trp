import os, socket
from trp import __path__
from os.path import join

VERBOSE = 0

DATADIR = os.path.join(__path__[0], 'data')

# archive (backed-up) results here
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')

# cache for temporary files
CACHEDIR = join(os.path.expanduser('~'), '.trp_cache')
if not os.path.exists(CACHEDIR): os.mkdir(CACHEDIR)

for l in [DATADIR, RESULTSDIR, CACHEDIR]:
    if not os.path.exists(l):
        if VERBOSE:
            print(f"Making {l}")
        os.mkdir(l)
