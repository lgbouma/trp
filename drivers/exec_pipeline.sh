#!/bin/bash

# Unpack the envvironment (with the trp package), and activate it
tar -xzf py38.tar.gz
export PYTHONPATH=$PWD/py38

# Run the main Python script; called as run_trp.py <sample_id>.
python3 run_trp.py $1

# Move and tarball output (log  files, pkl files, anything downloaded from
# MAST) to top level directory for condor to return.  Everything will be tar'd.
save_log=true
save_pkl=false
save_mastDownload=false

mkdir $_CONDOR_SCRATCH_DIR/joboutput_$1

if [ "$save_log" = true ]; then
    mv /srv/.trp_cache/rotperiod_finding/*$1*/*log $_CONDOR_SCRATCH_DIR/joboutput_$1/.
fi
if [ "$save_pkl" = true ]; then
    mv /srv/.trp_cache/rotperiod_finding/*$1*/*pkl $_CONDOR_SCRATCH_DIR/joboutput_$1/.
fi
if [ "$move_mastDownload" = true ]; then
    mv /srv/.trp_cache/rotperiod_finding/*$1*/mastDownload $_CONDOR_SCRATCH_DIR/joboutput_$1/.
fi

tar czf joboutput_$1.tar.gz $_CONDOR_SCRATCH_DIR/joboutput_$1
