# pre_circure

The workflow is used to map reads to contigs to generate a .paf files used as input for the `circure` pipeline.  

## What it does 
- Minimap2: maps long-read sequencing data to assembled contigs using high-quality long-read presets, producing PAF alignment files
- PAF post-processing: trims Minimap2 output to the core alignment fields (first 11 columns) to generate a standardized, lighter weight PAF file for downstream analysis


## Conda environments
The pipeline uses Snakemake per-rule conda environments. The rule reference an environment by name for legacy/reuse, and exported to .yaml file found in the `workflows` directory..
Naming convention for exported envs:
`environment_<env-name>.yaml`
(e.g. `environment_minimap2.yaml`)
