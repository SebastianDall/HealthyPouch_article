# HealthyPouch_article	

## Microbiome Analysis
This repository contains the code used for generating the microbiome analysis. 

data/pipelines/analysis

´data´ contains

`data` contains: 
<ul>
  <li>`data/methylation_plot/`: Contains the data needed to generate the plots in `analysis/methylation_plot.R`:</li>
  <li>`data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt`: MetaPhlAn taxonomic profile generated from short-read metagenomic data. </li>
  <li>`data/gyr_arg.csv`: Contains the data needed to generate Table 2.</li>
  <li>`data/mag_qual.tsv`: Contains the MIMAG quality classifications for generated MAGs.</li>
  <li>`data/par_arg.csv`: Contains the data needed to generate Table 2.</li>
  <li>`data/participant_metadata.csv`: Contains the metadata for the included participants. </li>
</ul>


`analysis` contains: 
<ul>
  <li>`data/methylation_plot/`: Contains the data needed to generate the plots in `analysis/methylation_plot.R`:</li>
  <li>`data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt`: </li>
  <li>`data/gyr_arg.csv`: </li>
  <li>`data/mag_qual.tsv`: </li>
  <li>`data/par_arg.csv`: </li>
  <li>`data/methylation_plot/`: </li>
</ul>


 - `pipelines` contains: 
 - [Annotate plasmid args](plasmids/arg_pipelie/README.md)   </li>
 - [Validate contig circularity](pipelines/circure/README.md)   </li>
 - `data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt`: </li>
 - `data/gyr_arg.csv`: </li>
 - `data/mag_qual.tsv`: </li>
 - `data/par_arg.csv`: </li>
 - `data/methylation_plot/`: </li>




Microbiome Analysis
This repository contains the code used for generating the microbiome analysis. All scripts for code generation can be found in src folder.

An external library is needed, which can be installed as:

remotes::install_github("SebastianDall/mplibrary")

Scripts are:

src/heatmaps.Rmd: Produces heatmaps of the relative abundance of the top most abundant species.
src/alpha-beta_diversity.Rmd: Produces alpha and beta diversity plots.
src/ordinationplot.Rmd: Produces ordination plots.
Richness was defined as species with a relative abundance >0 and alpha diversity was calculated using the Shannon diversity index. Patient sample similarity to donors were calculated using both Sørensen coefficient on relative abundances and Bray-Curtis on Hellinger transformed relative abundances. The similarity would be measured as the median similarity to donor samples received and also as median similarity to a donor sample from each donor not received. For the placebo group, similarity was calculated as the median similarity to a donor sample from each donor used in the FMT group.

## Prerequisites
All data was analyzed using R (4.1.0) and RStudio.

## Data
Sequencing data is available at PRJEB66493 & PRJEB80556. To rerun the workflow install snakemake and run the snakemake pipeline as:


## Cite 