# circure

This pipeline uses the circure workflow to validate circularity of contigs.

## What it does
- Uses .paf alignment files produced by the `pre_circure` pipeline as input
- Applies circure scripts to assess and validate contig circularity based on read-mapping patterns

## Running the pipeline
The pipeline is executed in two parts, indicated as `PART 1` and `PART 2` in the Snakefile.
Only one rule all should be active at a time.

## Conda environments 
Conda environments are provided as YAML files in the envs/ directory.
