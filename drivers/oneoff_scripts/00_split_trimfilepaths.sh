#!/bin/bash

# Set the input file name
input_file="trimfilepaths.csv"

# Set the number of lines per output file
lines_per_file=10000

# Create the output directory if it doesn't exist
output_dir="trimfilepaths_chunks"
mkdir -p "$output_dir"

# Split the input file into multiple files
split -l "$lines_per_file" -d -a 5 "$input_file" "$output_dir/chunk_"

# Rename the output files to have a .csv extension
for file in "$output_dir"/chunk_*; do
    mv "$file" "$file.csv"
done

echo "File splitting completed."
