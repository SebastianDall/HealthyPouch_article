library(tidyverse)
library(ampvis2)


# Function
filter_tax_species <- function(dataset) {
  taxonomic_levels <- strsplit(dataset$clade_name, "\\|")
  max_level <- max(sapply(taxonomic_levels, length))
  
  keep_rows <- logical(nrow(dataset))
  keep_rows[1] <- TRUE
  
  for (i in 2:nrow(dataset)) {
    levels <- taxonomic_levels[[i]]
    has_strain <- any(grepl("^t__", levels))
    
    if (has_strain) {
      keep_rows[i] <- TRUE
    }
  }
  
  filtered_dataset <- dataset[keep_rows, ]
  return(filtered_dataset)
}

# Load data
metadata <- read_delim("data/participant_metadata.csv") %>%  
  select(sample_barcode, id, health_status) 
metaphlan <- read_delim("data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt", delim = "\t", skip=1, show_col_types = FALSE) %>%  
  rename_with(~ str_remove(.x, "_NonHuman_Combined_Subsampled_2500000")) %>%  
  filter_tax_species() %>% 
  rename_at(vars(-1), ~ sub("_.*$", "", .))


## Prepare amp object
metaphlan_df <- metaphlan %>%  
  filter(clade_name == "UNCLASSIFIED" | str_detect(clade_name, "s_{2}")) %>%  
  mutate(OTU = paste0("OTU", row_number())) %>%  
  relocate(OTU) %>%  
  mutate(clade_name = if_else(clade_name == "UNCLASSIFIED", "k__UNCLASSIFIED|p__UNCLASSIFIED|c__UNCLASSIFIED|o__UNCLASSIFIED|f__UNCLASSIFIED|g__UNCLASSIFIED|s__UNCLASSIFIED|t__UNCLASSIFIED", clade_name)) %>%
  separate(clade_name, sep = "\\|", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                                             "Species", "Strain"))

otutable <- metaphlan_df %>% 
  select(-Kingdom:-Strain)

taxtable <- metaphlan_df %>% 
  slice(-1) %>% 
  relocate(OTU, .before="Kingdom") %>% 
  mutate(Species = Species %>%
           sub("^s__", "", .) %>%   
           gsub("_", " ", .)) 
  
amp_object <- amp_load(otutable = otutable, metadata = metadata, taxonomy = taxtable)


## Plot
rda <- amp_object %>%  
  amp_ordinate(
    type = "rda", 
    transform = "hellinger", 
    filter_species = 0,
    constrain = "health_status", 
    sample_color_by = "health_status",
    species_plot = T,
    species_nlabels = 5,
    species_label_taxonomy = "Species", 
    sample_point_size = 2) + 
  scale_color_manual(name = NULL, 
                     labels = c("SickPouch", "HealthyPouch", "UC", "NormalGut"),
                     values = c(  "#803E75", "#A6BDD7", "#CEA262",  "#817066"),
                     breaks = c("sick_pouch_incl", "healthy_pouch", "UC", "normal_gut")) + 
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"), 
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
        axis.ticks.x = element_line(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"), 
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8), 
        axis.title.x =  element_text(color = "black", size = 8), 
        axis.title.y =  element_text(color = "black", size = 8, margin = margin(r = 1)), 
        legend.position = "bottom",
        legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 10, , margin = margin(l = 0)), 
        legend.key.width  = unit(6, "pt"),
        legend.box.margin = margin(t = -15, r = 0, b = 0, l = 0),
        plot.title = element_text(face = "bold", size = 10),
        text = element_text(family = "Times New Roman"), 
        aspect.ratio = 1) +
  labs(title = "D) Redundancy Analysis") 

for (i in seq_along(rda$layers)) {
  g <- class(rda$layers[[i]]$geom)[1]
  if (g %in% c("GeomText", "GeomTextRepel")) {
    rda$layers[[i]]$aes_params$family <- "Times New Roman"
    rda$layers[[i]]$aes_params$size   <- 2.5}}
