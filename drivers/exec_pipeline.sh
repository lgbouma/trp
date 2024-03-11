#!/bin/bash

# Unpack the envvironment (with the trp package), and activate it
tar -xzf py38.tar.gz
export PYTHONPATH=$PWD/py38

# Run the main Python script; called as run_trp.py <sample_id>.
python3 run_trp.py $1

# Move and tarball output (log and pkl files, plus the MAST download) to top
# level directory for condor return.
mkdir $_CONDOR_SCRATCH_DIR/joboutput_$1
mv /srv/.trp_cache/rotperiod_finding/*$1*/*log $_CONDOR_SCRATCH_DIR/joboutput_$1/.
mv /srv/.trp_cache/rotperiod_finding/*$1*/*pkl $_CONDOR_SCRATCH_DIR/joboutput_$1/.
mv /srv/.trp_cache/rotperiod_finding/*$1*/mastDownload $_CONDOR_SCRATCH_DIR/joboutput_$1/.
tar czf joboutput_$1.tar.gz $_CONDOR_SCRATCH_DIR/joboutput_$1
