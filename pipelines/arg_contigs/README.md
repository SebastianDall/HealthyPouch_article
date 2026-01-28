# arg_contigs

The workflow is uses CARD's RGI main tool to predict ARGs within all contigs of a given sample. 

`RGI` github repo can be found [here](https://github.com/arpcard/rgi) and the `RGI main` help page [here](https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst)

## What it does
- Scans **all contigs** in each sample assembly using **CARD RGI (`rgi main`)** to predict ARGs.
- Automatically loads the **local CARD database JSON** (configured via `DATABASE_JSON`) using `rgi load --local`.
- Runs RGI with **DIAMOND alignments** (`-a DIAMOND`) and flags suited for metagenomic/contig data (`--low_quality`, `--include_nudge`, `--clean`).
- Produces per-sample outputs in `{OUTDIR}/{sample}/`:
  - `{sample}.txt` (tabular RGI results)
  - `{sample}.json` (JSON results; currently written by the workflow)

## Database download
The newest version of the CARD database was downloaded into `/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/databases/CARD_3.3.0` using the following command:
```bash
wget https://card.mcmaster.ca/latest/data
tar -xvf data
```

## Conda environments
The pipeline uses Snakemake per-rule conda environments. The rule reference an environment by name for legacy/reuse, and exported to .yaml file found in the `workflows` directory.
Naming convention for exported envs:
`environment_<env-name>.yaml`
(e.g. `environment_rgi.yaml`)
