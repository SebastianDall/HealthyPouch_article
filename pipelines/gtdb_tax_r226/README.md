# gtdb_tax_r226

Pipeline used to taxonomically classify bins, generated using the `pipelines/mmlong2-lite` pipeline.

## What it does
- Runs GTDB-Tk (v2.4.1) classify_wf to assign GTDB taxonomy to assemblies (MAGs/bins).

- Processes genomes in batch folders (batch_*) under the configured input directory, using a configurable file extension (e.g. fa, fasta).
    For each batch, produces GTDB-Tk summary tables for:
        Bacteria (bac120 marker set): gtdbtk.bac120.summary.tsv
        Archaea (ar53 marker set): gtdbtk.ar53.summary.tsv

- If no archaeal genomes are detected in a batch, writes a placeholder message to gtdbtk.ar53.summary.tsv to keep downstream steps consistent.


## Database download
The newest version of the GTDB (226) was downloaded into `path/to/database` using the following command:

```bash
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar -xvzf gtdbtk_data.tar.gz
```

## Conda environments
The pipeline uses Snakemake per-rule conda environments. The rule reference an environment by name for legacy/reuse, and exported to .yaml file found in the `workflows` directory.
Naming convention for exported envs:
`environment_<env-name>.yaml`
(e.g. `environment_gtdbtk_2.4.1.yaml`)
