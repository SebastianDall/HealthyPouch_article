# gtdb_tax

Pipeline used to taxonomically classify the single-sample-binning bins. 

### Location 
/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article

### Installation and conda environment generation
Newest version (2.4.0) of GTDB-tk was used.

Build conda env gtdbtk: 
```bash
mamba create -p /home/projects/cu_00014/people/albmol/_environments/gtdbtk --channel conda-forge --channel bioconda --channel defaults gtdbtk=2.4.0
ln -s /home/projects/cu_00014/people/albmol/_environments/gtdbtk ~/.conda/envs/
```

### Database download
The newest version of the GTDB (220) was downloaded into `/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/databases/gtdb_r220` using the following command:

```bash
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
```

### Pre-pipeline modifications
Prior to running the pipeline symbolic links to all bins were created using the following script:
```bash
scripts/create_symlinks_gtdb_tax.sh
```


### Pipeline generation 
To run the snakefile, the working directory has to be: 
```bash
/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article
```

The command for performing a dry run is as follows: 
```bash
snakemake -s _pipelines/gtdb_tax/workflows/Snakefile --configfile _pipelines/gtdb_tax/config/config.yaml -n 
```

For executing the snakefile, use the following command: 
```bash
snakemake -s _pipelines/gtdb_tax/workflows/Snakefile --configfile _pipelines/gtdb_tax/config/config.yaml --cores 40 --use-conda
```

### Execution script 
File path: 
```bash
../scripts/gtdb_tax_execute.sh
```


### Submitting to qsub
File path: 
```bash
qsub scripts/gtdb_tax_execute.sh
```
