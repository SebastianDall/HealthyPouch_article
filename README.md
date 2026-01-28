# HealthyPouch_article	

## Microbiome Analysis
This repository contains the code used for generating the microbiome analysis. 

`data/` contains: 
 - `data/methylation_plot/`: Contains the data needed to generate the plots in `analysis/methylation_plot.R`.
 - `data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt`: MetaPhlAn taxonomic profile generated from short-read metagenomic data. </li>
 - `data/gyr_arg.csv`: Contains the data needed to generate Table 2.</li>
 - `data/mag_qual.tsv`: Contains the MIMAG quality classifications for generated MAGs.</li>
 - `data/par_arg.csv`: Contains the data needed to generate Table 2.</li>
 - `data/participant_metadata.csv`: Contains the metadata for the included participants. </li>

`analysis` contains: 
 - `analysis/alpha_diversity.R`: script used to generate plot A and B in Figure 1.
 - `analysis/heatmap.R`: script used to generate plot C in Figure 1.
 - `analysis/methylation_plot.R`: script used to generate plot B in Figure 3 and Supplementary Figured 1 and 2.
 - `analysis/plasmid_plot.R`: script used to generate plot A in Figure 3 and Supplementary Figures 1 and 2.
 - `analysis/rda.R`: script used to generate plot D in Figure 1.

`pipelines/` contains workflows to: 
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

## Cite 
Please cite our preprint: 
