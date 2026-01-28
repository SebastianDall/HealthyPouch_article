# mmlong2-lite

Pipeline for mmlong2-lite workflow. Note that the assembly was previously generated using mmlong2-lite, and is hence skipped in this pipeline. 

## What it does
- Bins previously generated assemblies using the mmlong2-lite workflow

### Donwloading mmlong2-lite files
There are several files and folders on the mmlong2-lite github page: https://github.com/Serka-M/mmlong2-lite

```bash
wget -P pipelines/mmlong2-lite/ https://raw.githubusercontent.com/Serka-M/mmlong2-lite/main/src/mmlong2-lite
wget -P pipelines/mmlong2-lite/ https://raw.githubusercontent.com/Serka-M/mmlong2-lite/main/src/mmlong2-lite-config.yaml
wget -P pipelines/mmlong2-lite/ https://raw.githubusercontent.com/Serka-M/mmlong2-lite/main/src/mmlong2-lite.smk
```

### Conda environments 

The pipeline uses Snakemake per-rule conda environments. Some rules reference YAMLs via `config.yaml`, while others reference an environment by name for legacy/reuse.

All environment definitions (including exports for name-based envs) are collected in:
`pipelines/mmlong2-lite/envs`

Naming convention for exported envs:
`environment_<env-name>.yaml`
(e.g. `environment_checkm_genome_1.2.3.yaml`)

Rules referencing YAMLs via `config.yaml`, were downloaded from mmlong2-lite github as follows: 
```bash
wget -P pipelines/mmlong2-lite/envs/ https://raw.githubusercontent.com/Serka-M/mmlong2-lite/main/src/envs/env_1.yaml
wget -P pipelines/mmlong2-lite/envs/ https://raw.githubusercontent.com/Serka-M/mmlong2-lite/main/src/envs/env_2.yaml
wget -P pipelines/mmlong2-lite/envs/ https://raw.githubusercontent.com/Serka-M/mmlong2-lite/main/src/envs/env_6.yaml
```
