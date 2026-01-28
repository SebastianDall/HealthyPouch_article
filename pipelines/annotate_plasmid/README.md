# annotate_plasmid

Workflow used to annotate plasmids. 

## What it does
- **Prodigal**: predicts ORFs from input contigs
- **DIAMOND**: searches predicted proteins against **MobileOG-db**

## Prior to running pipeline
### Database download and modifications
Database
Prior to running the pipeline, mobile-OG was downloaded and converted to diamond file.

mobile-OG (prerelase 2.0.1-90[pre-release]) is downloaded into `path/to/database/mobileOGdb_1.6`:
```bash

wget "https://mobileogdb-downloads.s3.us-east-2.amazonaws.com/data-version-files/beatrix-1-6_v1_all.zip?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA3AL4C2VS6QHGRAFU%2F20250623%2Fus-east-2%2Fs3%2Faws4_request&X-Amz-Date=20250623T104132Z&X-Amz-Expires=900&X-Amz-SignedHeaders=host&X-Amz-Signature=49396da5a645b82c20cb4e70fff7747979e94415f190d40bbc7aea4eccb36723" \
     -O beatrix-1-6_v1_all.zip
unzip beatrix-1-6_v1_all.zip
```

Convert fasta version of database to diamond version:
```bash
cd path/to/database/mobileOGdb_1.6
mamba activate diamond
diamond makedb --in mobileOGdb_90-2.0.fasta -d mobileOGdb_90-2.0.dmnd
```

### Modify script
Change permissions of script before running: 
```bash
chmod +x pipelines/annotate_plasmid/workflows/mobileOGs-pl-kyanite.sh
```

## Installation and conda environment generation
Build conda env containing all necessary dependencies. 

```bash
mamba create -p path/to/environments/mobileOG-db --channel conda-forge --channel bioconda --channel defaults python=3.6.15
ln -s path/to/environments/mobileOG-db ~/.conda/envs/

mamba activate mobileOG-db
mamba install -c conda-forge biopython
mamba install -c bioconda prodigal
mamba install -c bioconda diamond
mamba install -c anaconda pandas
```

The pipeline uses Snakemake per-rule conda environments. The rule reference an environment by name for legacy/reuse, and exported to .yaml file found in the `workflows` directory..
Naming convention for exported envs:
`environment_<env-name>.yaml`
(e.g. `environment_mobileOG-db.yaml`)


