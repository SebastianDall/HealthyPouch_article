#!/bin/bash

# Define the directories
source_directory="/home/projects/cu_00014/data/sepseq_WGS/analysis/AMR_analysis/Nanopore_WGS/data/20_flye_assembly"
metadata_file="/home/projects/cu_00014/data/sepseq_WGS/analysis/AMR_analysis/amr-pred-exp/data/df_ab_filtered_wide.tsv"
destination_directory="data/sepseq_assemblies"

# Create destination directory if it doesn't exist
mkdir -p "$destination_directory"

# Read the sample_id column from the metadata file
# Assuming the sample_id column is the first column (change if needed)
sample_ids=$(awk -F'\t' 'NR > 1 {print $1}' "$metadata_file")

# Loop through each sample ID
for sample_id in $sample_ids; do
    # Define the source assembly file path
    assembly_file="$source_directory/$sample_id/assembly.fasta"

    # Check if the assembly file exists
    if [[ -f "$assembly_file" ]]; then
        # Create a subdirectory in the destination directory
        subdirectory="$destination_directory/$sample_id"
        mkdir -p "$subdirectory"

        # Create a symbolic link to the assembly file
        ln -s "$assembly_file" "$subdirectory/assembly.fasta"
        echo "Created symlink for $sample_id"
    else
        echo "Assembly file for $sample_id does not exist: $assembly_file"
    fi
done
