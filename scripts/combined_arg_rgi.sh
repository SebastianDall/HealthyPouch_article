#!/bin/bash

# Define the base directory containing the folders
BASE_DIR="/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/output/arg_rgi"
# Define the output directory for the combined DataFrame
OUTPUT_DIR="/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/output/arg_rgi_combined"
# Define the output file
OUTPUT_FILE="$OUTPUT_DIR/combined_arg_rgi.txt"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Initialize a flag to check if the header has been added
HEADER_ADDED=false

# Loop through each folder in the base directory
for folder in "$BASE_DIR"/*; do
    if [ -d "$folder" ]; then  # Check if it's a directory
        # Get the name of the folder to use as sample_barcode
        sample_barcode=$(basename "$folder")
        # Find the .txt file in the current folder
        txt_file="$folder/$sample_barcode.txt"
        
        # Check if the .txt file exists
        if [ -f "$txt_file" ]; then
            # If header has not been added, add it and set the flag
            if [ "$HEADER_ADDED" = false ]; then
                # Add sample_barcode header to the output file
                echo -e "sample_barcode\t$(head -n 1 "$txt_file")" > "$OUTPUT_FILE"
                HEADER_ADDED=true
            fi
            # Append data from the .txt file, excluding the header (first line)
            tail -n +2 "$txt_file" | awk -v barcode="$sample_barcode" '{print barcode "\t" $0}' >> "$OUTPUT_FILE"
        fi
    fi
done

echo "Combined data saved to $OUTPUT_FILE"

