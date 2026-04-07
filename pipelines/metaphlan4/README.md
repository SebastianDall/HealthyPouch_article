# MetaPhlAn4 pipeline

This pipeline will remove human reads by mapping to the HG38 genome and subset all samples to equal read depth before making a
MetaPhlAn4 taxonomy table. Everything can be run from the `Snakefile`, where it is important to set the path to the HG38 genome
database. The `SUBSET` parameter will determine how many reads the samples are subset to.

A tsv file is necessary to guide the pipeline. In the file it is specified where the R1 and R2 fastq file is for which all the
subsequent steps will be done. The pipeline can handle a sample having multiple illumina runs, which will concatenate the read files
before subsampling. These are simply extra entries in the file.


The pipeline will do the following:
- Quality trim and filter reads using `fastp`
- Remove read mapping to the human genome using `bowtie2`
- Concatenate quality-trimmed, non-human read files if multiple runs have been done.
- Subsample concatenated read files to a specific depth using `seqtk`
- Create the MetaPhlAn4 taxonomy profile


