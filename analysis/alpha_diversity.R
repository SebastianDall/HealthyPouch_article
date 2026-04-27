library(tidyverse)
library(ampvis2)
library(vegan)
library(ggpubr)


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
    select(sample_barcode, health_status)
readcount <- read_delim("./data/read_count.txt", col_names = c("sample_barcode", "library_id", "num_reads"))
# metaphlan <- read_delim("./data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt", delim = "\t", skip=1, show_col_types = FALSE) %>%  
#     rename_with(~ str_remove(.x, "_NonHuman_Combined_Subsampled_full")) %>%  
#     mutate(clade_name = if_else(clade_name == "UNCLASSIFIED", "k__UNCLASSIFIED|p__UNCLASSIFIED|c__UNCLASSIFIED|o__UNCLASSIFIED|f__UNCLASSIFIED|g__UNCLASSIFIED|s__UNCLASSIFIED", clade_name)) %>%
#     separate_wider_delim(clade_name, delim="|", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), too_few = "align_start") %>%
#     filter(!is.na(Species)) %>%
#     filter(is.na(Strain)) 
metaphlan <- read_delim("data/MetaPhlAn_4.1.0_Combined_NonHuman_Subsampled_full_profile.txt", delim = "\t", skip=1, show_col_types = FALSE) %>%  
    rename_with(~ str_remove(.x, "_NonHuman_Combined_Subsampled_full")) %>%  
    mutate(clade_name = if_else(clade_name == "UNCLASSIFIED", "k__UNCLASSIFIED|p__UNCLASSIFIED|c__UNCLASSIFIED|o__UNCLASSIFIED|f__UNCLASSIFIED|g__UNCLASSIFIED|s__UNCLASSIFIED", clade_name)) %>%
    separate_wider_delim(clade_name, delim="|", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), too_few = "align_start") %>%
    filter(!is.na(Species)) %>%
    filter(is.na(Strain)) 
# filter_tax_species() %>% 
# rename_at(vars(-1), ~ sub("_.*$", "", .))


readcount_sum <- readcount %>%
    select(sample_barcode, num_reads) %>%
    group_by(sample_barcode) %>%
    summarise(num_reads = sum(num_reads))

sample_barcode_filtered <- readcount_sum %>%
    filter(num_reads >= 2500000) %>% pull(sample_barcode)

 
## Dataframe
readcount_plot_df <- readcount_sum %>%
    filter(sample_barcode %in% sample_barcode_filtered) %>%
    right_join(metadata, by = "sample_barcode") %>%
    mutate(health_status = factor(health_status, levels = c("sick_pouch_incl", "healthy_pouch", "UC", "normal_gut")))

## Dataframe
readcount_plot_removed_df <- readcount_sum %>%
    filter(!sample_barcode %in% sample_barcode_filtered) %>%
    right_join(metadata, by = "sample_barcode") %>%
    mutate(health_status = factor(health_status, levels = c("sick_pouch_incl", "healthy_pouch", "UC", "normal_gut")))


## Statistics
comparisons <- tribble(
  ~group_code, ~g1,               ~g2,
  "sp_hp",     "sick_pouch_incl",  "healthy_pouch",
  "sp_uc",     "sick_pouch_incl",  "UC",
  "sp_ng",     "sick_pouch_incl",  "normal_gut",
  "hp_uc",     "healthy_pouch",    "UC",
  "hp_ng",     "healthy_pouch",    "normal_gut",
  "uc_ng",     "UC",              "normal_gut"
)

## Statistics
wilcox_pair_reads <- function(df, g1, g2) {
    d <- df %>% filter(health_status %in% c(g1, g2))
    tibble(
           g1 = g1,
           g2 = g2,
           group = str_c(g1, " ~ ", g2),
           sig_test_p_value = wilcox.test(num_reads ~ health_status, data = d)$p.value
           ) %>%
    mutate(sig_dif = if_else(sig_test_p_value <= 0.05, "yes", "no"))
}

readcount_stats <- comparisons %>%
    mutate(res = map2(g1, g2, ~ wilcox_pair_reads(readcount_plot_df, .x, .y))) %>%
    select(group_code, res) %>%
    unnest(res) %>%
    mutate(
           p_val_adj = p.adjust(sig_test_p_value, method = "bonferroni"),
           p_val_adj_sim = case_when(
                                     p_val_adj < 0.0001 ~ "p<0.0001",
                                     p_val_adj < 0.001  ~ "p<0.001",
                                     p_val_adj < 0.01   ~ "p<0.01",
                                     p_val_adj < 0.05   ~ "p<0.05",
                                     p_val_adj > 0.05   ~ "p>0.05",
                                     TRUE               ~ NA_character_
           )
    )

## Plot
readcount_plot <- ggplot(readcount_plot_df, aes(x = health_status, y = num_reads)) +
    geom_point(aes(color = health_status), position = position_jitter(width = 0.09), alpha = 0.7, size = 2) +
    geom_hline(yintercept = 2500000, linetype = "dashed") +
    geom_point(data = readcount_plot_removed_df, color = "black", position = position_jitter(width = 0.09), alpha = 0.7, size = 2) +
    scale_x_discrete(labels = c("SickPouch", "HealthyPouch", "UC", "NormalGut")) +
    scale_color_manual(name = NULL,
                       values = c("#803E75", "#A6BDD7", "#CEA262", "#817066")) +
    theme(panel.background = element_rect(fill = "grey97"),
          panel.grid.major = element_line(color = "grey85"),
          panel.grid.minor = element_line(color = "grey85", linetype = "dotted"),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.line = element_line(color = "black", linewidth = 0.1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          legend.position = "none",
          axis.text.x = element_text(color = "black", size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(size = 8),
          text = element_text(family = "Times New Roman"),
          plot.title = element_text(face = "bold", size = 10),
          aspect.ratio = 1) +
    labs(x = "", y = "Read Count", title = "Read Count")

readcount_plot


# Richness plot
metaphlan <- metaphlan %>% select(Kingdom:Strain, sample_barcode_filtered)
view(metaphlan)

## Dataframe
richness_df <- metaphlan %>%  
    filter(Species != "s__UNCLASSIFIED") %>%
    select(-c(Kingdom:Genus, Strain)) %>%
    column_to_rownames(var = "Species") %>% 
    t() %>%  
    specnumber() %>%  
    enframe() %>%  
    rename(sample_barcode = name, 
           richness = value) %>% 
    right_join(metadata, by = "sample_barcode") %>%  
    mutate(health_status = factor(health_status, levels = c("sick_pouch_incl", "healthy_pouch", "UC", "normal_gut")))


wilcox_pair <- function(df, g1, g2) {
  d <- df %>% filter(health_status %in% c(g1, g2))
  tibble(
    g1 = g1,
    g2 = g2,
    group = str_c(g1, " ~ ", g2),
    sig_test_p_value = wilcox.test(richness ~ health_status, data = d)$p.value
  ) %>%
    mutate(sig_dif = if_else(sig_test_p_value <= 0.05, "yes", "no"))
}

richness_stats <- comparisons %>%
  mutate(res = map2(g1, g2, ~ wilcox_pair(richness_df, .x, .y))) %>%
  select(group_code, res) %>%
  unnest(res) %>%
  mutate(
    p_val_adj = p.adjust(sig_test_p_value, method = "bonferroni"),
    p_val_adj_sim = case_when(
      p_val_adj < 0.0001 ~ "p<0.0001",
      p_val_adj < 0.001  ~ "p<0.001",
      p_val_adj < 0.01   ~ "p<0.01",
      p_val_adj < 0.05   ~ "p<0.05",
      p_val_adj > 0.05   ~ "p>0.05",
      TRUE               ~ NA_character_
    )
  )


## Plot
richness_plot <- ggplot(richness_df, aes(x = health_status, y = richness)) +
  geom_point(aes(color = health_status), position = position_jitter(width = 0.09), alpha = 0.7, size = 2) + 
  scale_x_discrete(labels = c("SickPouch", "HealthyPouch", "UC", "NormalGut")) + 
  scale_color_manual(name = NULL, 
                     values = c(  "#803E75", "#A6BDD7", "#CEA262",  "#817066")) + 
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"), 
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
        axis.ticks.x = element_line(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "none", 
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8), 
        axis.title.y = element_text(size = 8), 
        text = element_text(family = "Times New Roman"), 
        plot.title = element_text(face = "bold", size = 10), 
        aspect.ratio = 1) + 
  labs(x = "", y = "Richness") + 
  geom_bracket(xmin = "sick_pouch_incl", xmax = "healthy_pouch", y.position = 175,
               label = "p < 0.05", label.size = 8/.pt, family = "Times New Roman") +
  geom_bracket(xmin = "healthy_pouch", xmax = "UC", y.position = 275,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") + 
  geom_bracket(xmin = "UC", xmax = "normal_gut", y.position = 310,
               label = "p < 0.05", label.size = 8/.pt, family = "Times New Roman") + 
  geom_bracket(xmin = "sick_pouch_incl", xmax = "UC", y.position = 240,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") + 
  geom_bracket(xmin = "healthy_pouch", xmax = "normal_gut", y.position = 345,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") + 
  geom_bracket(xmin = "sick_pouch_incl", xmax = "normal_gut", y.position = 380,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") +
  labs(title = "A) Richness")


# Shannon diversity

## Dataframe
metaphlan_df <- metaphlan %>%  
  mutate(OTU = paste0("OTU", row_number())) %>%  
  relocate(OTU)

otutable <- metaphlan_df %>% 
  select(-Kingdom:-Strain)

taxtable <- metaphlan_df %>% 
  # slice(-1) %>% 
  relocate(OTU, .before="Kingdom")

amp_object <- amp_load(otutable = otutable, metadata = metadata, taxonomy = taxtable)

shannon_div_df <- amp_alpha_diversity(amp_object, measure = "shannon", richness = FALSE, rarefy = NULL) %>% 
  mutate(health_status = factor(health_status, levels = c("sick_pouch_incl", "healthy_pouch", "UC", "normal_gut")))


## Statistics
wilcox_pair <- function(df, g1, g2) {
  d <- df %>% filter(health_status %in% c(g1, g2))
  tibble(
    g1 = g1,
    g2 = g2,
    group = str_c(g1, " ~ ", g2),
    sig_test_p_value = wilcox.test(Shannon ~ health_status, data = d)$p.value
  ) %>%
    mutate(sig_dif = if_else(sig_test_p_value <= 0.05, "yes", "no"))
}

shannon_div_stats <- comparisons %>%
  mutate(res = map2(g1, g2, ~ wilcox_pair(shannon_div_df, .x, .y))) %>%
  select(group_code, res) %>%
  unnest(res) %>%
  mutate(
    p_val_adj = p.adjust(sig_test_p_value, method = "bonferroni"),
    p_val_adj_sim = case_when(
      p_val_adj < 0.0001 ~ "p<0.0001",
      p_val_adj < 0.001  ~ "p<0.001",
      p_val_adj < 0.01   ~ "p<0.01",
      p_val_adj < 0.05   ~ "p<0.05",
      p_val_adj > 0.05   ~ "p>0.05",
      TRUE               ~ NA_character_
    )
  )


## Plot
shannon_div_plot <- ggplot(shannon_div_df, aes(x = health_status, y = Shannon)) +
  geom_point(aes(color = health_status), position = position_jitter(width = 0.09), alpha = 0.7, size = 2) + scale_x_discrete(labels = c("SickPouch", "HealthyPouch", "UC", "NormalGut")) + 
  scale_color_manual(name = NULL, 
                     values = c(  "#803E75", "#A6BDD7", "#CEA262",  "#817066")) + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(panel.background = element_rect(fill="grey97"),
        panel.grid.major = element_line(color = "grey85"), 
        panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
        axis.ticks.x = element_line(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "none", 
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8), 
        axis.title.y = element_text(size = 8, margin = margin(r = 20)), 
        text = element_text(family = "Times New Roman"), 
        plot.title = element_text(face = "bold", size = 10), 
        aspect.ratio = 1) + 
  labs(x = "", y = "Shannon Diversity Index") +
  geom_bracket(xmin = "healthy_pouch", xmax = "UC", y.position = 4.8,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") + 
  geom_bracket(xmin = "sick_pouch_incl", xmax = "UC", y.position = 5.7,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") + 
  geom_bracket(xmin = "healthy_pouch", xmax = "normal_gut", y.position = 6.3,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") + 
  geom_bracket(xmin = "sick_pouch_incl", xmax = "normal_gut", y.position = 6.6,
               label = "p < 0.0001", label.size = 8/.pt, family = "Times New Roman") +
  labs(title = "B) Shannon Diversity")

shannon_div_plot
