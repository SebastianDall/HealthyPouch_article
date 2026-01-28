# plasmid_iden

This pipeline runs four different plasmid tools, on contigs. 

## What it does
Runs a multi-tool contig classification workflow per sample to identify plasmid signals from an input FASTA.

- Runs multiple tools for plasmid classification, specificall: 
    ViralVerify (1.1)
    Plasmer (23.04.20)
    PlasmidHunter
    geNomad

## Database download

***viralverify:***
Downloaded into `path/to/databases/viralverify_hmms` using the following command:

```bash
curl -L -J -O "https://figshare.com/ndownloader/files/17904323?private_link=f897d463b31a35ad7bf0"
gunzip nbc_hmms.hmm.gz
```

***plasmer:***
Downloaded into `path/to/databases/plasmer` using the following command:
```bash
curl -L -J -O "https://zenodo.org/records/7030675/files/customizedKraken2DB.tar.xz?download=1"
tar -xf customizedKraken2DB.tar.xz
``` 

```bash
curl -L -J -O "https://zenodo.org/records/7030675/files/plasmerMainDB.tar.xz?download=1"
tar -xf plasmerMainDB.tar.xz
```

***PlasmidHunter***
Part of environment:
```bash
curl -L -J -O "https://zenodo.org/records/7030675/files/plasmerMainDB.tar.xz?download=1"
tar -xf plasmerMainDB.tar.xz
```

***genomad***
Downloaded into `path/to/databases/genomad` using as advsed by genomad [github](https://github.com/apcamargo/genomad).


## Conda environments
The pipeline uses Snakemake per-rule conda environments. The rule reference an environment by name for legacy/reuse, and exported to .yaml file found in the `workflows` directory.
Naming convention for exported envs:
`environment_<env-name>.yaml`
(e.g. `environment_viralverify.yaml`)





