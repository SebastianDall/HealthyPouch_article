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

ra_df <- amp_heatmap(amp_object,
                               group_by = "sample_barcode", 
                               tax_aggregate = "Species",
                               plot_values = FALSE,
                               normalise = FALSE, 
                               textmap = TRUE,
                               tax_show = 15) 


## Aggregate and cluster dataframe
aggregated_rab <- ra_df %>%
  rownames_to_column(var = "species") %>% 
  pivot_longer(cols = ASYM00000001MP:UCPr00000032MP, names_to = "sample_barcode", values_to = "ra") %>% 
  left_join(metadata, by = "sample_barcode") %>% 
  group_by(health_status, species) %>%
  summarise(mean_abundance = mean(ra, na.rm = TRUE)) %>%
  ungroup() %>%  
  mutate(
    facet_label = case_when(
      health_status == "sick_pouch_incl" ~ "SickPouch",
      health_status == "healthy_pouch"   ~ "HealthyPouch",
      health_status == "UC"              ~ "UC",
      health_status == "normal_gut"      ~ "NormalGut",
      TRUE                               ~ NA_character_)) %>%  
  mutate(facet_label = factor(facet_label, levels = c("SickPouch", "HealthyPouch", "UC", "NormalGut")))

matrix <- aggregated_rab %>% 
  pivot_wider(
    names_from = species,           
    values_from = mean_abundance) %>%  
  select(-facet_label) %>%  
  column_to_rownames(var = "health_status")

species_dist <- dist(t(matrix), method = "euclidean")  
species_hclust <- hclust(species_dist, method = "ward.D2")

species_order <- species_hclust$order
species_labels <- species_hclust$labels  

aggregated_rab_ordered <- aggregated_rab %>%
  mutate(
    species = str_replace_all(species, "_", " "),
    species = factor(species, levels = str_replace_all(species_labels[species_order], "_", " "))
  ) %>%  
  mutate(mean_abundance = round(mean_abundance, digits=1))


## Plot
heatmap <- ggplot(aggregated_rab_ordered, aes(x = health_status, y = species, fill = mean_abundance)) +
  geom_tile() +
  geom_text(aes(label = round(mean_abundance, 2)), size = 8/.pt, color = "black", family = "Times New Roman") + 
  facet_grid(. ~ facet_label, scales = "free_x", space = "free") + 
  labs(y = "Species", fill = "% relative\n abundance", x = "") +
  scale_y_discrete(limits = levels(aggregated_rab_ordered$species), drop = TRUE) + 
  scale_fill_gradient(low = "gray90", high = "#8E7A47") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, face = "italic", color="black"),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(), 
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5), 
    panel.grid = element_line(color = "grey80"),
    panel.background = element_rect(fill = "grey97"),
    plot.background = element_rect(fill = "transparent", color = "transparent"), 
    strip.text = element_text(size = 8), 
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = "bottom", 
    legend.key.height = unit(5, "mm"),
    legend.key.width  = unit(5, "mm"),
    legend.box.margin = margin(t = -15, r = 0, b = 0, l = 0), 
    plot.title = element_text(face = "bold", size = 10), 
    text = element_text(family = "Times New Roman")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(title = "C) Relative Abundance")
heatmap
