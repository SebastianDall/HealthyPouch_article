library(tidyverse)
library(ampvis2)
library(vegan)


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
pca <- amp_object %>%  
  amp_ordinate(
    type = "pca", 
    # transform = "hellinger", 
    filter_species = 0,
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
  labs(title = "PCA") 

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


ggpubr::ggarrange(pca, rda, nrows = 1)

ggsave(filename = "./figures/pca-rda_short_read.png", device = "png", dpi = "retina", bg = "white")

for (i in seq_along(rda$layers)) {
  g <- class(rda$layers[[i]]$geom)[1]
  if (g %in% c("GeomText", "GeomTextRepel")) {
    rda$layers[[i]]$aes_params$family <- "Times New Roman"
    rda$layers[[i]]$aes_params$size   <- 2.5}}




# ── 1. Prepare data for vegan ──────────────────────────────────────
#p Extract OTU matrix (samples × species) and Hellinger-transform
otu_mat <- t(amp_object$abund)                 # samples as rows
otu_hell <- decostand(otu_mat, method = "hellinger")

view(otu_hell)

# Group factor
health <- amp_object$metadata$health_status

# ── 2. RDA permutation test (anova.cca) ────────────────────────────
# This directly tests whether constrained axes explain more than chance
rda_model <- rda(otu_hell ~ health_status, data = amp_object$metadata)

# Global test (all axes)
set.seed(42)
rda_global <- anova.cca(rda_model, permutations = 9999)
print(rda_global)

# Per-axis test
set.seed(42)
rda_byaxis <- anova.cca(rda_model, by = "axis", permutations = 50)
print(rda_byaxis)

# Variance explained by the constraint
R2 <- RsquareAdj(rda_model)
cat("R² =", round(R2$r.squared, 4),
    "\nAdj. R² =", round(R2$adj.r.squared, 4), "\n")

# ── 3. PERMANOVA (adonis2) ─────────────────────────────────────────
# Gold-standard test for multivariate group differences
set.seed(42)
permanova_euc <- adonis2(otu_hell ~ health_status,
                         data = amp_object$metadata,
                         method = "euclidean",   # Euclidean on Hellinger = Hellinger distance
                         permutations = 999)
print(permanova_euc)

# Also on Bray-Curtis (raw abundances) for robustness
set.seed(42)
permanova_bray <- adonis2(otu_mat ~ health_status,
                          data = amp_object$metadata,
                          method = "bray",
                          permutations = 999)
print(permanova_bray)

# Pairwis PERMANOVA:
groups <- levels(factor(health))
pairs <- combn(groups, 2, simplify = FALSE)

pairwise_results <- do.call(rbind, lapply(pairs, function(pair) {
                                              idx <- health %in% pair
                                              set.seed(42)
                                              res <- adonis2(otu_hell[idx, ] ~ health_status,
                                                             data = amp_object$metadata[idx, ],
                                                             method = "euclidean", permutations = 999)
                                              data.frame(
                                                         Group1 = pair[1], Group2 = pair[2],
                                                         F_stat = res$F[1], R2 = res$R2[1], p_value = res$`Pr(>F)`[1]
                                              )
                                             }))

pairwise_results$p_adj <- p.adjust(pairwise_results$p_value, method = "BH")
print(pairwise_results)

# ── 4. Betadisper – homogeneity of dispersions ────────────────────
# PERMANOVA assumes equal dispersion; this checks that assumption
dist_hell <- dist(otu_hell, method = "euclidean")
betadisp <- betadisper(dist_hell, health)

set.seed(42)
betadisp_test <- permutest(betadisp, permutations = 999)
print(betadisp_test)

# Pairwise dispersion comparisons
print(TukeyHSD(betadisp))

tukey <- TukeyHSD(betadisp)

tukey_d <- tukey$group %>%
    as_tibble() %>%
    mutate(group = row.names(tukey$group)) %>%
    rename(p_adj = `p adj`) %>%
    pivot_longer(cols = diff:upr, names_to = "boundary", values_to = "value")

tukey_d %>%
ggplot(aes(x = value, y = group)) +
geom_line() +
geom_point(shape = "|", size = 5) +
geom_label(data = tukey_d %>% filter(boundary == "diff"), aes(label = round(p_adj, 2)), nudge_y = 0.3) +
geom_vline(xintercept = 0, linetype = "dashed") +
labs(title ="TukeyHSD on betadisp test")

# ── 5. Pairwise PERMANOVA ─────────────────────────────────────────
# Which groups actually differ from each other?
groups <- levels(factor(health))
pairs <- combn(groups, 2, simplify = FALSE)

pairwise_results <- do.call(rbind, lapply(pairs, function(pair) {
                                              idx <- health %in% pair
                                              set.seed(42)
                                              res <- adonis2(otu_hell[idx, ] ~ health_status,
                                                             data = amp_object$metadata[idx, ],
                                                             method = "euclidean", permutations = 999)
                                              data.frame(
                                                         Group1 = pair[1], Group2 = pair[2],
                                                         F_stat = res$F[1], R2 = res$R2[1], p_value = res$`Pr(>F)`[1]
                                              )
                                             }))

# FDR correction
pairwise_results$p_adj <- p.adjust(pairwise_results$p_value, method = "BH")
print(pairwise_results)

# ── 6. ANOSIM (supplementary) ─────────────────────────────────────
set.seed(42)
anosim_res <- anosim(dist_hell, health, permutations = 999)
print(anosim_res)

