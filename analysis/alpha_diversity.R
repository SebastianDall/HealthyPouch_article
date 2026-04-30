library(tidyverse)
library(ampvis2)
library(vegan)
library(ggpubr)


# Load data
metadata <- read_delim("data/participant_metadata.csv") %>%  
    select(sample_barcode, health_status)
#readcount <- read_delim("./data/read_count.txt", col_names = c("sample_barcode", "library_id", "num_reads"))
metaphlan <- read_delim("./data/MetaPhlAn_4.1.0_NonHuman_Subsampled_2500000_profile.txt", delim = "\t", skip=1, show_col_types = FALSE) %>%  
    rename_with(~ str_remove(.x, "_NonHuman_Combined_Subsampled_2500000")) %>%  
    mutate(clade_name = if_else(clade_name == "UNCLASSIFIED", "k__UNCLASSIFIED|p__UNCLASSIFIED|c__UNCLASSIFIED|o__UNCLASSIFIED|f__UNCLASSIFIED|g__UNCLASSIFIED|s__UNCLASSIFIED", clade_name)) %>%
    separate_wider_delim(clade_name, delim="|", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), too_few = "align_start") %>%
    filter(!is.na(Species)) %>%
    filter(is.na(Strain)) 


# readcount_sum <- readcount %>%
#     select(sample_barcode, num_reads) %>%
#     group_by(sample_barcode) %>%
#     summarise(num_reads = sum(num_reads))


# sample_barcode_filtered <- readcount_sum %>%
#     filter(num_reads >= 2500000) %>% pull(sample_barcode)

# ## Dataframe
# readcount_plot_df <- readcount_sum %>%
#     filter(sample_barcode %in% sample_barcode_filtered) %>%
#     right_join(metadata, by = "sample_barcode") %>%
#     mutate(health_status = factor(health_status, levels = c("sick_pouch_incl", "healthy_pouch", "UC", "normal_gut")))

# ## Dataframe
# readcount_plot_removed_df <- readcount_sum %>%
#     filter(!sample_barcode %in% sample_barcode_filtered) %>%
#     right_join(metadata, by = "sample_barcode") %>%
#     mutate(health_status = factor(health_status, levels = c("sick_pouch_incl", "healthy_pouch", "UC", "normal_gut")))


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
# 
# ## Statistics
# wilcox_pair_reads <- function(df, g1, g2) {
#     d <- df %>% filter(health_status %in% c(g1, g2))
#     tibble(
#            g1 = g1,
#            g2 = g2,
#            group = str_c(g1, " ~ ", g2),
#            sig_test_p_value = wilcox.test(num_reads ~ health_status, data = d)$p.value
#            ) %>%
#     mutate(sig_dif = if_else(sig_test_p_value <= 0.05, "yes", "no"))
# }

# readcount_stats <- comparisons %>%
#     mutate(res = map2(g1, g2, ~ wilcox_pair_reads(readcount_plot_df, .x, .y))) %>%
#     select(group_code, res) %>%
#     unnest(res) %>%
#     mutate(
#            p_val_adj = p.adjust(sig_test_p_value, method = "bonferroni"),
#            p_val_adj_sim = case_when(
#                                      p_val_adj < 0.0001 ~ "p<0.0001",
#                                      p_val_adj < 0.001  ~ "p<0.001",
#                                      p_val_adj < 0.01   ~ "p<0.01",
#                                      p_val_adj < 0.05   ~ "p<0.05",
#                                      p_val_adj > 0.05   ~ "p>0.05",
#                                      TRUE               ~ NA_character_
#            )
#     )

# ## Plot
# readcount_plot <- ggplot(readcount_plot_df, aes(x = health_status, y = num_reads)) +
#     geom_point(aes(color = health_status), position = position_jitter(width = 0.09), alpha = 0.7, size = 2) +
#     geom_hline(yintercept = 2500000, linetype = "dashed") +
#     geom_point(data = readcount_plot_removed_df, color = "black", position = position_jitter(width = 0.09), alpha = 0.7, size = 2) +
#     scale_x_discrete(labels = c("SickPouch", "HealthyPouch", "UC", "NormalGut")) +
#     scale_color_manual(name = NULL,
#                        "values = c("#803E75", "#A6BDD7", "#CEA262", "#817066")
#                        values = c(  "#803E75", "#4f6985", "#CEA262",  "#817066") +
#     theme(panel.background = element_rect(fill = "grey97"),
#           panel.grid.major = element_line(color = "grey85"),
#           panel.grid.minor = element_line(color = "grey85", linetype = "dotted"),
#           axis.ticks.x = element_line(color = "black"),
#           axis.ticks.y = element_line(color = "black"),
#           axis.line = element_line(color = "black", linewidth = 0.1),
#           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
#           plot.background = element_rect(fill = "transparent", color = "transparent"),
#           legend.position = "none",
#           axis.text.x = element_text(color = "black", size = 8),
#           axis.text.y = element_text(color = "black", size = 8),
#           axis.title.y = element_text(size = 8),
#           text = element_text(family = "Times New Roman"),
#           plot.title = element_text(face = "bold", size = 10),
#           aspect.ratio = 1) +
#     labs(x = "", y = "Read Count", title = "Read Count")
# 
# readcount_plot


# Richness plot
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

add_sig_bracket_richness <- function(xmin, xmax, y.position) {
  
  p_label <- richness_stats %>%
    filter(
      (g1 == xmin & g2 == xmax) |
        (g1 == xmax & g2 == xmin)
    ) %>%
    pull(p_val_adj_sim) %>%
    str_replace("p", "p[adj] ")
  
  if (length(p_label) == 0 || p_label == "p[adj] >0.05") {
    return(NULL)
  }
  
  geom_bracket(
    xmin = xmin,
    xmax = xmax,
    y.position = y.position,
    label = p_label,
    label.size = 8/.pt,
    family = "Times New Roman",
    type = "expression"
  )
}


## Plot
richness_plot <- ggplot(richness_df, aes(x = health_status, y = richness)) +
  geom_point(aes(color = health_status), position = position_jitter(width = 0.09), alpha = 0.7, size = 2) + 
  scale_x_discrete(labels = c("SickPouch", "HealthyPouch", "UC", "NormalGut")) + 
  scale_color_manual(name = NULL, 
                     #values = c(  "#803E75", "#A6BDD7", "#CEA262",  "#817066")
                     values = c("#803E75", "#4f6985", "#CEA262",  "#817066")) + 
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
  ylim(0,360)+
  labs(x = "", y = "Richness") +
  add_sig_bracket_richness("sick_pouch_incl", "healthy_pouch", 175) +
  add_sig_bracket_richness("healthy_pouch", "UC", 275) +
  add_sig_bracket_richness("sick_pouch_incl", "UC", 240) +
  add_sig_bracket_richness("healthy_pouch", "normal_gut", 310) +
  add_sig_bracket_richness("sick_pouch_incl", "normal_gut", 343) +
  add_sig_bracket_richness("UC", "normal_gut", 343) + 
  labs(title = "A) Richness")

richness_plot


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


add_sig_bracket_shannon <- function(xmin, xmax, y.position) {
  
  p_label <- shannon_div_stats %>%
    filter(
      (g1 == xmin & g2 == xmax) |
        (g1 == xmax & g2 == xmin)
    ) %>%
    pull(p_val_adj_sim) %>%
    str_replace("p", "p[adj] ")
  
  if (length(p_label) == 0 || p_label == "p[adj] >0.05") {
    return(NULL)
  }
  
  geom_bracket(
    xmin = xmin,
    xmax = xmax,
    y.position = y.position,
    label = p_label,
    label.size = 8/.pt,
    family = "Times New Roman",
    type = "expression"
  )
}


## Plot
shannon_div_plot <- ggplot(shannon_div_df, aes(x = health_status, y = Shannon)) +
  geom_point(aes(color = health_status), position = position_jitter(width = 0.09), alpha = 0.7, size = 2) + scale_x_discrete(labels = c("SickPouch", "HealthyPouch", "UC", "NormalGut")) + 
  scale_color_manual(name = NULL, 
                     #values = c(  "#803E75", "#A6BDD7", "#CEA262",  "#817066")
                     values = c(  "#803E75", "#4f6985", "#CEA262",  "#817066")) + 
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
  ylim(0,7.2) + 
  labs(x = "", y = "Shannon Diversity Index") +
  add_sig_bracket_shannon("healthy_pouch", "UC", 4.8) + 
  add_sig_bracket_shannon("sick_pouch_incl", "UC", 5.7) + 
  add_sig_bracket_shannon("healthy_pouch", "normal_gut", 6.3) + 
  add_sig_bracket_shannon("sick_pouch_incl", "normal_gut", 6.8) +
  add_sig_bracket_shannon("UC", "normal_gut", 5.3) +
  add_sig_bracket_shannon("sick_pouch_incl", "healthy_pouch", 2) + 
  labs(title = "B) Shannon Diversity")
  
shannon_div_plot

