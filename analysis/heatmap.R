library(tidyverse)
library(ampvis2)
library(Maaslin2)

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

metaphlan <- read_delim("./data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt", delim = "\t", skip=1, show_col_types = FALSE) %>%  
    rename_with(~ str_remove(.x, "_NonHuman_Combined_Subsampled_full")) %>%  
    mutate(clade_name = if_else(clade_name == "UNCLASSIFIED", "k__UNCLASSIFIED|p__UNCLASSIFIED|c__UNCLASSIFIED|o__UNCLASSIFIED|f__UNCLASSIFIED|g__UNCLASSIFIED|s__UNCLASSIFIED", clade_name)) %>%
    separate_wider_delim(clade_name, delim="|", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), too_few = "align_start") %>%
    filter(!is.na(Species)) %>%
    filter(is.na(Strain)) 
# metaphlan <- read_delim("data/MetaPhlAn_4.1.0_Combined_NonHuman_Subsampled_full_profile.txt", delim = "\t", skip=1, show_col_types = FALSE) %>%  
#     rename_with(~ str_remove(.x, "_NonHuman_Combined_Subsampled_full")) %>%  
#     mutate(clade_name = if_else(clade_name == "UNCLASSIFIED", "k__UNCLASSIFIED|p__UNCLASSIFIED|c__UNCLASSIFIED|o__UNCLASSIFIED|f__UNCLASSIFIED|g__UNCLASSIFIED|s__UNCLASSIFIED", clade_name)) %>%
#     separate_wider_delim(clade_name, delim="|", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), too_few = "align_start") %>%
#     filter(!is.na(Species)) %>%
#     filter(is.na(Strain)) 
#
#
# readcount <- read_delim("./data/read_count.txt", col_names = c("sample_barcode", "library_id", "num_reads"))
# readcount_sum <- readcount %>%
#     select(sample_barcode, num_reads) %>%
#     group_by(sample_barcode) %>%
#     summarise(num_reads = sum(num_reads))
#
# sample_barcode_filtered <- readcount_sum %>%
#     filter(num_reads >= 2500000) %>% pull(sample_barcode)
#
# metaphlan <- metaphlan %>% select(Kingdom:Strain, sample_barcode_filtered)


## Prepare amp object
metaphlan_df <- metaphlan %>%  
  mutate(OTU = paste0("OTU", row_number())) %>%  
  relocate(OTU)

otutable <- metaphlan_df %>% 
  select(-Kingdom:-Strain)

otutable %>% select(-OTU) %>% summarise(across(everything(), ~sum(.)))

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


# Prepare input: species × samples matrix (species as rows, samples as columns)
# Your otutable already has this structure, just need to transpose
otu_for_maaslin <- amp_object$abund %>% t() %>% as.data.frame()

# Metadata
meta_for_maaslin <- amp_object$metadata %>%
    select(-sample_barcode) %>%
    mutate(health_status = factor(health_status,
                                  levels = c("normal_gut", "healthy_pouch", "sick_pouch_incl", "UC")))

# Run MaAsLin2 — reference is normal_gut
fit <- Maaslin2(
                input_data = otu_for_maaslin / 100,
                input_metadata = meta_for_maaslin,
                output = "results/maaslin2",
                fixed_effects = "health_status",
                reference = c("health_status,healthy_pouch"),
                normalization = "CLR",        # already relative abundance, but MaAsLin2 renormalizes
                transform = "NONE",
                min_abundance = 0.01,         # filter very rare taxa
                min_prevalence = 0.1,         # present in at least 10% of samples
                analysis_method = "LM",
                plot_heatmap = FALSE,
                plot_scatter = FALSE
)

# Significant results
sig_results <- fit$results %>%
    filter(qval < 0.05) %>%
    left_join(taxtable %>% select(OTU:Strain) %>% rename(feature = OTU)) %>%
    arrange(qval)



sig_results %>% filter(str_detect(Species, "catenibacter"))
summary(sig_results)

sig_results %>% filter(row_number() < 10)

top15_features <- sig_results %>%
  filter(qval < 0.05) %>%
  group_by(feature) %>%
  summarise(min_qval = min(qval), .groups = "drop") %>%
  arrange(min_qval) %>%
  head(15) %>%
  pull(feature)

# Prepare plot data
plot_data <- sig_results %>%
  filter(feature %in% top15_features) %>%
  mutate(
    log10q = -log10(qval),
    comparison = case_when(
      value == "healthy_pouch"   ~ "HealthyPouch",
      value == "normal_gut"   ~ "NormalGut",
      value == "sick_pouch_incl" ~ "SickPouch",
      value == "UC"              ~ "UC"
    ),
    sig_label = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01  ~ "**",
      qval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )

# Order species by mean absolute coefficient for cleaner visual
species_order <- plot_data %>%
  group_by(Species) %>%
  summarise(mean_abs_coef = mean(abs(coef)), .groups = "drop") %>%
  arrange(mean_abs_coef) %>%
  pull(Species)

plot_data <- plot_data %>%
  mutate(
    Species = factor(Species, levels = species_order),
    comparison = factor(comparison, levels = c("SickPouch", "UC", "NormalGut"))
  )

# Bubble plot
da_plot <- ggplot(plot_data, aes(x = comparison, y = Species)) +
  geom_point(aes(size = log10q, fill = coef), shape = 21, stroke = 0.3, color = "grey30") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "Coefficient\n(log fold change)"
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(q))),
    range = c(2, 10),
    breaks = c(2, 5, 10, 20)
  ) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.75, family = "Times New Roman") +
  labs(
    x = "Group (vs NormalGut)",
    y = NULL,
    title = "Top 15 differentially abundant species"
  ) +
  theme(
    panel.background = element_rect(fill = "grey97"),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8, face = "italic"),
    axis.title.x = element_text(color = "black", size = 8),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.1),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    legend.text = element_text(color = "black", size = 8),
    legend.title = element_text(color = "black", size = 8),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 10),
    text = element_text(family = "Times New Roman")
  )

da_plot

ggsave("figures/differential_abundance_top15.png", da_plot,
       device = "png", dpi = "retina", bg = "white",
       width = 7, height = 5)


# ── MaAsLin2: HealthyPouch vs SickPouch ───────────────────────────
meta_for_maaslin_pouch <- amp_object$metadata %>%
  select(-sample_barcode) %>%
  mutate(health_status = factor(health_status,
                                levels = c("healthy_pouch", "sick_pouch_incl", "normal_gut", "UC")))

fit_pouch <- Maaslin2(
  input_data = otu_for_maaslin / 100,
  input_metadata = meta_for_maaslin_pouch,
  output = "results/maaslin2_pouch_ref",
  fixed_effects = "health_status",
  reference = c("health_status,healthy_pouch"),
  normalization = "CLR",
  transform = "NONE",
  min_abundance = 0.01,
  min_prevalence = 0.1,
  analysis_method = "LM",
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

sig_results_pouch <- fit_pouch$results %>%
  filter(value == "sick_pouch_incl", qval < 0.05) %>%
  left_join(taxtable %>% select(OTU:Strain) %>% rename(feature = OTU)) %>%
  arrange(qval)

top15_pouch <- sig_results_pouch %>%
  arrange(qval) %>%
  head(15) %>%
  pull(feature)

plot_data_pouch <- sig_results_pouch %>%
  filter(feature %in% top15_pouch) %>%
  mutate(
    log10q = -log10(qval),
    sig_label = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01  ~ "**",
      qval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )

species_order_pouch <- plot_data_pouch %>%
  arrange(coef) %>%
  pull(Species)

plot_data_pouch <- plot_data_pouch %>%
  mutate(Species = factor(Species, levels = species_order_pouch))

da_plot_pouch <- ggplot(plot_data_pouch, aes(x = "SickPouch vs HealthyPouch", y = Species)) +
  geom_point(aes(size = log10q, fill = coef), shape = 21, stroke = 0.3, color = "grey30") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "Coefficient\n(log fold change)"
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(q))),
    range = c(2, 10),
    breaks = c(2, 5, 10, 20)
  ) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.75, family = "Times New Roman") +
  labs(
    x = NULL,
    y = NULL,
    title = "Top 15 differentially abundant species\n(SickPouch vs HealthyPouch)"
  ) +
  theme(
    panel.background = element_rect(fill = "grey97"),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8, face = "italic"),
    axis.title.x = element_text(color = "black", size = 8),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.1),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    legend.text = element_text(color = "black", size = 8),
    legend.title = element_text(color = "black", size = 8),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 10),
    text = element_text(family = "Times New Roman")
  )

da_plot_pouch

ggsave("figures/differential_abundance_sick_vs_healthy_pouch.png", da_plot_pouch,
       device = "png", dpi = "retina", bg = "white",
       width = 6, height = 5)


# ── Heatmap with MaAsLin2 significance markers ───────────────────
sig_for_heatmap <- fit$results %>%
  left_join(taxtable %>% select(OTU, Species) %>% rename(feature = OTU)) %>%
  mutate(
    species = str_replace_all(Species, "_", " "),
    facet_label = case_when(
      value == "healthy_pouch"   ~ "HealthyPouch",
      value == "normal_gut"   ~ "NormalGut",
      value == "sick_pouch_incl" ~ "SickPouch",
      value == "UC"              ~ "UC"
    ),
    sig_label = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01  ~ "**",
      qval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  ) %>%
  filter(!is.na(facet_label)) %>%
  select(species, facet_label, qval, sig_label)

heatmap_sig_data <- aggregated_rab_ordered %>%
  left_join(sig_for_heatmap, by = c("species", "facet_label")) %>%
  mutate(
    sig_label = replace_na(sig_label, ""),
    tile_label = if_else(sig_label == "",
                         as.character(round(mean_abundance, 2)),
                         paste0(round(mean_abundance, 2), sig_label)),
    facet_label = factor(facet_label, levels = c("SickPouch", "HealthyPouch", "UC", "NormalGut"))
  )

heatmap_sig <- ggplot(heatmap_sig_data, aes(x = health_status, y = species, fill = mean_abundance)) +
  geom_tile() +
  geom_text(aes(label = tile_label), size = 8/.pt, color = "black", family = "Times New Roman") +
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
  labs(title = "C) Relative Abundance",
       caption = "Significance vs HealthyPouch (MaAsLin2 CLR): * q<0.05, ** q<0.01, *** q<0.001")
heatmap_sig

ggsave("figures/heatmap_with_significance.png", heatmap_sig,
       device = "png", dpi = "retina", bg = "white",
       width = 7, height = 5)
