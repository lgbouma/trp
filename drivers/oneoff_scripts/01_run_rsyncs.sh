#!/bin/bash

# Set the directory containing the chunk files
chunk_dir="trimfilepaths_chunks"

# Set the remote host and username
remote_host="luke@wh1.caltech.edu"

# Set the log directory
log_dir="logs_v2"

# Create the log directory if it doesn't exist
mkdir -p "$log_dir"

# Loop through the chunk files
for chunk_file in "$chunk_dir"/chunk_*.csv; do

    # Extract the chunk number from the file name
    chunk_number=$(echo "$chunk_file" | sed -E 's/.*chunk_([0-9]+)\.csv/\1/')
    
    # Set the log file name
    log_file="$log_dir/rsynclog_$chunk_number.log"
    
    # Run the rsync command in the foreground with the chunk file and log the output
    rsync -chavzP --files-from="$chunk_file" "$remote_host":/ . &> "$log_file"
    
    echo "Started rsync for $chunk_file. Logging to $log_file."
done

echo "All rsync commands have been run."
