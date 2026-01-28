# virsorter2

## What it does
- Runs a virsorter2 on contigs to identify viral sequences.

## Conda environments
The pipeline uses Snakemake per-rule conda environments. The rule reference an environment by name for legacy/reuse, and exported to .yaml file found in the `workflows` directory.
Naming convention for exported envs:
`environment_<env-name>.yaml`
(e.g. `environment_virsorter2.yaml.yaml`)

