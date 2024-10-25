#!/bin/bash

# Define source and target directories
SRC_DIR="/home/projects/cu_00014/data/healthypouch/binning/20240224_single_sample_binning"
DEST_DIR="/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/20240224_single_sample_binning_bins"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through all subdirectories (samples) in the source directory
for sample_dir in "$SRC_DIR"/*; do
    # Check if it's a directory
    if [ -d "$sample_dir" ]; then
        sample=$(basename "$sample_dir")
        bins_dir="$sample_dir/results/bins"

        # Check if the bins directory exists
        if [ -d "$bins_dir" ]; then
            echo "Processing $sample: bins directory found."
            
            # Loop through the files in the bins directory
            for bin_file in "$bins_dir"/*; do
                if [ -f "$bin_file" ]; then
                    # Create a symbolic link in the destination directory
                    ln -s "$bin_file" "$DEST_DIR/$(basename "$bin_file")"
                    echo "Link created for $bin_file"
                fi
            done
        else
            echo "No bins directory found for $sample."
        fi
    fi
done
