import os, socket
from trp import __path__
from os.path import join

DATADIR = os.path.join(__path__[0], 'data')

# archive (backed-up) results here
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')

for l in [DATADIR, RESULTSDIR]:
    if not os.path.exists(l):
        print(f"Making {l}")
        os.mkdir(l)
