#!/bin/bash

# Purpose: bulk rsync <500pc QLP light curves from wh1 to 
# openscigrid architecture.
#
# Usage:
# 1) From $DATA/QLP, run rsync_500pc_make_dirs.py (via singularity container)
# 2) From $DATA/QLP, ./rsync_500pc_from_wh1.sh &> get500pc_v0.log &

# Function to perform bulk rsync based on the CSV file
rsync_bulk() {

    csvpath="s0001_to_s0055_QLP_ticid_path_merged_parallax_gt_2mas.csv"

    Ntot=$(tail -n +2 "$csvpath" | wc -l)

    # Boolean option to control directory creation
    create_dirs=false

    # Create a temporary file to store the list of remote file paths
    tmp_file=$(mktemp)

    # Extract the remote file paths from the CSV and write them to the temporary file
    tail -n +2 "$csvpath" | cut -d',' -f2 | sed 's|^|luke@wh1.caltech.edu:/ar1/TESS/QLP/|' > "$tmp_file"

    # Create directories for all the files
    if [ "$create_dirs" = true ]; then
        while IFS=',' read -r ticid relpath; do
            if [[ $ticid != "ticid" ]]; then
                dirname=$(dirname "$relpath")
                if [[ ! -d $dirname ]]; then
                    mkdir -p "$dirname"
                fi
            fi
        done < "$csvpath"
    fi

    # Perform rsync using the --files-from option
    rsync -chavzP --no-motd --omit-dir-times --files-from="$tmp_file" luke@wh1.caltech.edu:/ar1/TESS/QLP/ ./

    # Remove the temporary file
    rm "$tmp_file"
}

# Run the rsync_bulk function
rsync_bulk
