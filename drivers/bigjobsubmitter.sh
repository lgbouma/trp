#!/bin/bash

# Purpose:
#    Queue jobs on OSG in a given volume slice, up to the 10k job max.
# Usage:
#   `./bigjobsubmitter.sh 0 100` (to submit trp chunks from 0 to 100pc)

# Check if the required arguments are provided
if [ $# -ne 2 ]; then
  echo "Usage: $0 <start_index> <end_index>"
  exit 1
fi

start_index=$1
end_index=$2

# Validate the start and end indices
if [ $start_index -lt 0 ] || [ $end_index -le $start_index ]; then
  echo "Invalid start or end index. Make sure start_index is non-negative and end_index is greater than start_index."
  exit 1
fi

# Get the generic username
username=$(whoami)

# Iterate over the specified range
for ((i=start_index; i<end_index; i+=5)); do
  start=$i
  end=$((i+5))
  volumeslice="${start}to${end}pc"
  
  # Check the number of jobs submitted by the user
  job_count=$(condor_q | grep "Total for $username" | awk '{print $4}')
  echo `condor_q | grep "Total for $username"`
  echo $job_count
  
  if [ -z "$job_count" ]; then
    job_count=0
  fi
  
  # Wait until the job count is below 8000 before submitting the next job
  while [ $job_count -ge 8000 ]; do
    echo "Waiting for job count to drop below 8000. Current count: $job_count"
    sleep 60  # Wait for 60 seconds before checking again
    job_count=$(condor_q | grep "Total for $username" | awk '{print $4}')
  done
  
  # Submit the condor job
  condor_submit "volumeslice=$volumeslice" condor_submit.sub
  echo "Submitted job for volumeslice: $volumeslice"
  sleep 2
done
