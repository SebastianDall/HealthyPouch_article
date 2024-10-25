SOURCE_BASE="/home/projects/cu_00014/data/healthypouch/binning/20240224_single_sample_binning"
TARGET_BASE="/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/20240224_single_sample_binning_assemblies"

# Create symbolic links for all <sample> directories
for sample_dir in "$SOURCE_BASE"/*; do
    # Check if it's a directory
    if [ -d "$sample_dir" ]; then
        # Extract the sample name
        sample=$(basename "$sample_dir")

        # Define the target directory for the symbolic link
        target_dir="$TARGET_BASE/$sample"
        
        # Create the target directory if it doesn't exist
        mkdir -p "$target_dir"

        # Define the source and target paths for the symbolic link
        source_file="$SOURCE_BASE/$sample/tmp/eukfilt/asm_pol_lenfilt_eukfilt.fasta"
        target_file="$target_dir/asm_pol_lenfilt_eukfilt.fasta"

        # Create the symbolic link if the source file exists
        if [ -f "$source_file" ]; then
            ln -sf "$source_file" "$target_file"
            echo "Created symlink for $sample"
        else
            echo "Source file does not exist for $sample, skipping..."
        fi
    fi
done
