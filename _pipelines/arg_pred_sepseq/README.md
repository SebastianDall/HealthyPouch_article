# arg_pred_sepseq

This pipeline is a copy of the `arg_contigs` pipeline, modified to predict args in sepseq genomes (that have phenotypic metadata). 

The workflow is uses CARD's RGI main tool to predict ARGs within all contigs of a given sample. 

`RGI` github repo can be found [here](https://github.com/arpcard/rgi) and the `RGI main` help page [here](https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst)


### Location 
/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/_pipelines


### Pre-pipeline modifications --> OBS EVT SLET
Prior to running the pipeline symbolic links to assemblies were created using the following script: 
```bash
scripts/create_symlinks_arg_pred_sepseq.sh
```

### Pipeline generation 
To run the snakefile, the working directory has to be: 
```bash
/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article
```

The command for performing a dry run is as follows: 
```bash
snakemake -s _pipelines/arg_pred_sepseq/workflows/Snakefile --configfile _pipelines/arg_pred_sepseq/config/config.yaml -n 
```

For executing the snakefile, use the following command: 
```bash
snakemake -s _pipelines/arg_pred_sepseq/workflows/Snakefile --configfile _pipelines/arg_pred_sepseq/config/config.yaml --cores 40 --use-conda 
```


### Execution script 
File path: 
```bash
../scripts/arg_pred_sepseq_execute.sh
```


### Submitting to qsub
File path: 
```bash
qsub scripts/arg_pred_sepseq_execute.sh
```
