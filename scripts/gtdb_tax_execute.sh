#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=cu_00014 -A cu_00014
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N 2024-10-25_gtdb_tax
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e logs/gtdb_tax_241025
#PBS -o logs/gtdb_tax_241025
### Only send mail when job is aborted or terminates abnormally
#PBS -M albertehm@bio.aau.dk
#PBS -m abe
### Number of nodes
#PBS -l nodes=1:thinnode:ppn=40
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 12 hours)
#PBS -l walltime=2:00:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load tools
module load perl/5.20.2

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in separate script and execute from here
module load snakemake/8.20.5

mkdir -p logs/gtdb_tax_241025

snakemake -s _pipelines/gtdb_tax/workflows/Snakefile --use-conda --cores 40 --keep-going 

#snakemake --jobs 50 \
#    --snakefile _pipelines/rgi_bwt/workflows/Snakefile \
#    --use-conda --conda-frontend conda \
#    --latency-wait 60 \
#    --cluster "qsub -W group_list=cu_00014 -A cu_00014 -l walltime={params.walltime} -l nodes=1:{params.nodetype}:ppn={threads} -l mem={resources.mem} -m n -o logs/RGI-bwt/ -e logs/RGI-bwt/"

