# HealthyPouch_article	

## Microbiome Analysis
This repository contains the code used for generating the microbiome analysis. 

`data` contains: 
 - `data/methylation_plot/`: Contains the data needed to generate the plots in [`methylation_plot.R`](analysis/methylation_plot.R).
 - `data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt`: MetaPhlAn taxonomic profile generated from short-read metagenomic data. </li>
 - `data/gyr_arg.csv`: Contains the data needed to generate Table 2.</li>
 - `data/mag_qual.tsv`: Contains the MIMAG quality classifications for generated MAGs.</li>
 - `data/par_arg.csv`: Contains the data needed to generate Table 2.</li>
 - `data/participant_metadata.csv`: Contains the metadata for the included participants. </li>

`analysis` contains: 
 - [`alpha_diversity.R`](analysis/alpha_diversity.R): script used to generate plot A and B in Figure 1.
 - [`heatmap.R`](analysis/heatmap.R): script used to generate plot C in Figure 1.
 - [`mag_qual.R`](analysis/mag_qual.R): script used to make figure 2. This one requires several files that can be found in the gtdb release 226 folder. Paths in the script will have to be adjusted.
 - [`methylation_plot.R`](analysis/methylation_plot.R): script used to generate plot B in Figure 3 and Supplementary Figured 1 and 2.
 - [`plasmid_plot.R`](analysis/plasmid_plot.R): script used to generate plot A in Figure 3 and Supplementary Figures 1 and 2.
 - [`rda.R`](analysis/rda.R): script used to generate plot D in Figure 1.

`pipelines` contains workflows to: 
 - [MiMAG gene idenfication pipeline](pipelines/mimag_pipeline/README.md) 
 - [Assembly stats](pipelines/stats/README.md) 
 - [Annotate plasmids](pipelines/annotate_plasmid/README.md) 
 - [Identify ARGs](pipelines/arg_contigs/README.md)
 - [Validate contig circularity](pipelines/circure/README.md)
 - [Taxonomically classify MAGs](pipelines/gtdb_tax_r226/README.md)
 - [Bin assembled contigs](pipelines/mmlong2-lite/README.md)
 - [Identify plasmids](pipelines/plasmid_iden/README.md)
 - [Prepare for contig circularity validation](pipelines/pre_circure/README.md)
 - [Identify virus](pipelines/viral_iden/README.md)


## Prerequisites
All data was analyzed using R (4.1.0) and RStudio.

## Data
All DNA sequences have been deposited at the European Nucleotide Archive under the accession number PRJEB98069. Additionally, we incorporated data from previous studies, which are publicly available under the accession numbers PRJEB80556, PRJEB6649

A script can be found with the accession and overview of all uploaded data at: `data/filereport_read_run_PRJEB98069.tsv`.
The files `data/hp_included_participants.csv`, `data/PaPr_metadata.tsv`, `data/PaOP_metadata.tsv`, and `data/donor_metadata.tsv` contains a `sample_barcode` column that matches the `sample_title` of the download script.

## Cite 
Please cite our preprint: 
