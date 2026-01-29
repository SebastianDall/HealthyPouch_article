library(data.table)
library(tidyverse)
library(ggbeeswarm)

gtdb_me <- fread("./data/mefinder.tsv")
gtdb_me_transposons <- gtdb_me %>% 
  rename(closest_genome_reference = bin) %>% 
  # filter(type %in% c("composite transposon", "ime", "ice")) %>% 
  # filter(type == "composite transposon") %>% 
  group_by(closest_genome_reference) %>% 
  summarise(
    n_transposons = n()
  )

# Here clean just means the ssu_all_r226_ref.tsv where the delimeter has been set to tab.
gtdb_ssu <- fread("./data/ssu_all_r226_ref_clean.tsv", header = F, col.names = c("gtdb_genome_reference_id", "ssu_accession", "classification"), sep = "\t", sep2 = NA)


gtdb_bac_rep <- fread("./data/bac120_metadata_r226.tsv")
gtdb_arc_rep <- fread("./data/ar53_metadata_r226.tsv")

gtdb_rep <- bind_rows(gtdb_bac_rep, gtdb_arc_rep) %>% 
  filter(gtdb_representative == "t")

gtdb_redun <- bind_rows(gtdb_bac_rep, gtdb_arc_rep) %>% 
  group_by(gtdb_taxonomy) %>% 
  summarise(cluster_size = n()) %>% 
  mutate(
    is_singleton = ifelse(cluster_size == 1, T,F),
    # closest_genome_reference = str_remove(gtdb_genome_representative, "^[A-Z]+_")
  )

genes <- fread("./bakta/annotation_summary.tsv")
genes <- genes %>% 
  separate_wider_delim(bin, delim = ".", names = c("sample", "bin", "round", "bin_id"), cols_remove = F) %>% 
  mutate(
    sample = ifelse(sample == "NP-don-Asym00000007MP", "NP-Asym00000007MP", sample),
    bin = case_when(
      sample == "NP-Asym00000007MP" ~ str_replace(bin, "NP-don-Asym00000007MP", "NP-Asym00000007MP"),
      TRUE ~ bin,
    ),
    ssu_16s_len = pmax(bakta_16s_len, barrnap_16s_len),
    ssu_23s_len = pmax(bakta_23s_len, barrnap_23s_len),
    ssu_5s_len = pmax(bakta_5s_len, barrnap_5s_len)
  )




papr88 <- fread("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/mmlong2_lite/NP-PaPr00000088MP/tmp/binning/round_1/checkm2/quality_report.tsv") %>% 
  mutate(sample = "NP-PaPr00000088MP")

papr88_cov <- fread("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/mmlong2_lite/NP-PaPr00000088MP/tmp/binning/round_1/cov.tsv") %>% 
  rename(contig = contigname) %>% 
  left_join(fread("./files/NP-PaPr00000088MP_contig_bin.tsv", col.names = c("contig", "bin"))) %>% 
  filter(!is.na(bin)) %>% 
  group_by(bin) %>% 
  summarise(
    total_length = sum(contigLen),
    cov = sum(totalAvgDepth * contigLen) / total_length
  ) %>% 
  mutate(sample = "NP-PaPr00000088MP")


papr11 <- fread("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/mmlong2_lite/NP-PaPr00000011MP/tmp/binning/round_1/checkm2/quality_report.tsv") %>% 
  mutate(sample = "NP-PaPr00000011MP")

papr11_cov <- fread("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/mmlong2_lite/NP-PaPr00000011MP/tmp/binning/round_1/cov.tsv") %>% 
  rename(contig = contigname) %>% 
  left_join(fread("./files/NP-PaPr00000011MP.contig_bin.1.1.tsv", col.names = c("contig", "bin"))) %>% 
  filter(!is.na(bin)) %>% 
  group_by(bin) %>% 
  summarise(
    total_length = sum(contigLen),
    cov = sum(totalAvgDepth * contigLen) / total_length
  ) %>% 
  mutate(sample = "NP-PaPr00000011MP")

p <- bind_rows(papr88, papr11) %>% 
  rename(
    bin = Name,
    contamination_checkm2 = Contamination,
    completeness_checkm2 = Completeness,
    n_contigs = Total_Contigs,
    genome_size = Genome_Size,
    contig_n50 = Contig_N50
  ) %>% 
  select(sample, bin, contamination_checkm2, completeness_checkm2, n_contigs, genome_size, contig_n50) %>% 
  left_join(bind_rows(papr88_cov, papr11_cov)) %>% 
  select(-total_length) %>% 
  drop_na()
  
checkm <- tibble(
  sample = genes$sample
) %>% 
  distinct() %>% 
  filter(!sample %in% c("NP-PaPr00000088MP", "NP-PaPr00000011MP")) %>% 
  mutate(
    quality_report = map(.x = sample, ~fread(paste0("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/mmlong2_lite/",.x, "/results/", .x, "_bins.tsv")))
  )

checkm <- checkm %>% 
  unnest(quality_report) %>% 
  rename(n_contigs = contigs) %>% 
  bind_rows(p) %>% 
  mutate(
    is_circular = str_detect(bin, "bin\\.c")
  )

mimag <- checkm %>% 
  left_join(genes %>% select(-c(round, bin_id))) %>% 
  mutate(
    operon = case_when(
      (bakta_16s > 0 | barrnap_16s >0) & (bakta_23s > 0 | barrnap_23s > 0) & (bakta_5s > 0 | barrnap_5s > 0) ~ 1,
      TRUE ~ 0
    ),
    mimag = case_when(
      completeness_checkm2 > 90 & contamination_checkm2 < 5 & operon == 1 & custom_trna_uniq >= 18 ~ "HQ",
      completeness_checkm2 >= 50 & contamination_checkm2 < 10 ~ "MQ",
      TRUE ~ "LQ"
    )
  )


mimag %>% 
  group_by(mimag) %>% 
  summarise(
    n = n()
  )
write_delim(mimag, "data/mimag_bin_quality.tsv", delim ="\t")

# 1. Distinct bins: 12390 -> 12299
# 2. remove LQ bins: 12299-3=12296
# 3. Bin NP-PaPr00000011MP.bin.1.2 is in quality_report, but does not have a bin -1 = 12295


# New MIMAG stats
mimag_stats_distinct <- mimag %>% 
  distinct(bin, .keep_all = TRUE) %>% 
  filter(bin != "NP-PaPr00000011MP.bin.1.2")


mimag_stats_distinct  %>% 
  group_by(mimag) %>% 
  summarise(
    n = n()
  )



gtdb_mag_tax <- tibble(
  batch = list.dirs("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/output/mmlong2_lite_gtdb/", recursive = F, full.names = T)
) %>% 
  mutate(
    bac = map(.x = batch, ~(fread(file.path(.x, "gtdbtk.bac120.summary.tsv")))),
    arc = map(.x = batch, ~(fread(file.path(.x, "gtdbtk.ar53.summary.tsv")))),
    batch_id = basename(batch)
  )

gtdb_arc <- gtdb_mag_tax %>% 
  filter(map_dbl(.x = arc, ~(nrow(.x))) > 0) %>% 
  select(arc, batch_id) %>% 
  unnest(arc) %>%
  rename(bin = user_genome)

gtdb_bac <- gtdb_mag_tax %>% 
  filter(map_dbl(.x = bac, ~(nrow(.x))) > 0) %>% 
  select(bac, batch_id) %>% 
  unnest(bac) %>% 
  mutate(
    closest_genome_reference_radius = as.numeric(na_if(closest_genome_reference_radius, "N/A")),
    closest_genome_ani = as.numeric(na_if(closest_genome_ani, "N/A")),
    closest_genome_af = as.numeric(na_if(closest_genome_af, "N/A")),
    closest_placement_radius = as.numeric(na_if(closest_placement_radius, "N/A")),
    closest_placement_ani = as.numeric(na_if(closest_placement_ani, "N/A")),
    closest_placement_af = as.numeric(na_if(closest_placement_af, "N/A"))
  ) %>% 
  rename(bin = user_genome) %>% 
  relocate(batch_id)



gtdb <- bind_rows(gtdb_arc, gtdb_bac)

####
gtdb_redo <- gtdb %>%  
  filter(batch_id %in% c("batch_redo", "batch_redo_unrun")) %>%  
  mutate(sample_barcode = sub("\\..*$", "", bin)) %>% 
  relocate(sample_barcode)

gtdb_num <- gtdb %>%  
  filter(!(batch_id %in% c("batch_redo", "batch_redo_unrun"))) %>%
  mutate(sample_barcode = sub("\\..*$", "", bin)) %>%  
  filter(!(sample_barcode %in% gtdb_redo$sample_barcode))


gtdb_fuzzy <- fread("gtdb_novel/gtdbtk.ani_closest.tsv") %>% 
  rename(bin = user_genome) %>% 
  filter(satisfies_gtdb_circumscription_criteria == "True") %>% 
  select(bin, reference_genome, reference_taxonomy)


gtdb <- bind_rows(gtdb_redo, gtdb_num) %>% 
  mutate(bin = str_remove(bin, "don-")) %>% 
  left_join(gtdb_fuzzy) %>%
  mutate(
    closest_genome_reference = case_when(
      !is.na(reference_genome) ~ reference_genome,
      T ~ closest_genome_reference
    ),
    closest_genome_taxonomy = case_when(
      !is.na(reference_taxonomy) ~ reference_taxonomy,
      T ~ closest_genome_reference
    ),
    classification = case_when(
      !is.na(reference_taxonomy) ~ reference_taxonomy,
      T ~ classification
    ),
  ) %>%
  separate(classification, sep = ";", into = c("domain", "phylum", "class", "order", "family", "genus", "species"), remove = F) %>% 
  select(-c(reference_genome, reference_taxonomy))
####

load_file <- function(file, ...) {
  if (file.exists(file)) {
    fread(file, ...)
  } else {
    tibble()
  }
}


# me <- mimag_stats_distinct %>% select(bin) %>% 
#   mutate(
#     me = map(.x = bin, ~(load_file(file.path("bakta", .x, "mobile_element_finder", paste0(.x, "_mef.csv.csv.csv")), skip=5)))
#   )
# 
# me_f <- me %>% 
#   filter(map_dbl(me, nrow) > 0) %>% 
#   mutate(me = map(me, ~(.x %>% mutate(mge_no = as.character(mge_no))))) %>% 
#   unnest(me)
# 
# me_transposons <- me_f %>% 
#   # filter(type %in% c("composite transposon", "ime", "ice")) %>% 
#   group_by(bin) %>% 
#   summarise(
#     n_transposons = n()
#   )

gtdb_mimag_stats_distinct <- mimag_stats_distinct %>% 
  # left_join(me_transposons) %>% 
  left_join(gtdb) %>% 
  relocate(sample_barcode)
gtdb %>% filter(bin %in% c("NP-ASYM00000007MP.bin.1.10", "NP-Asym00000007MP.bin.1.10"))

x <- gtdb_mimag_stats_distinct %>% filter(is.na(classification))
#write_delim(gtdb_mimag_stats_distinct, file = "output/genome_stats.tsv", delim ="\t")

# gtdb_mimag_stats_distinct %>% 
#   select(closest_genome_reference) %>% 
#   drop_na() %>% 
#   distinct() %>% 
#   write_delim("data/study_ref.csv")

novelty <- gtdb_mimag_stats_distinct %>% 
  select(bin, classification:species, mimag) %>% 
  pivot_longer(-c(bin, mimag, classification), names_to = "taxonomy", values_to = "value") %>% 
  mutate(
    value_str = str_remove(value, "([a-z]__)"),
    value_int = ifelse(value_str == "", 1, 0)
  ) %>% 
  mutate(
    taxonomy = factor(taxonomy, levels = c("domain", "phylum", "class", "order", "family", "genus", "species"), labels = c("domain", "phylum", "class", "order", "family", "genus", "species"))
  )

novelty %>% 
  filter(!is.na(value_str)) %>% 
  group_by(taxonomy, mimag) %>% 
  summarise(
    dis = n_distinct(value)
  )%>% 
  print(n= 25)

novelty %>% 
  filter(mimag != "LQ") %>% 
  filter(!is.na(value_str)) %>% 
  group_by(taxonomy) %>% 
  summarise(
    dis = n_distinct(value)
  )%>% 
  print(n= 25)

novelty %>% 
  filter(mimag != "LQ") %>% 
  group_by(taxonomy, mimag) %>% 
  drop_na() %>%
  summarise(
    novel = sum(value_int),
    total = n()
  ) %>% 
  arrange(taxonomy) %>% 
  print(n= 25)


gtdb_rep_sel <- gtdb_rep %>%
  select(gtdb_type_designation_ncbi_taxa, genome_size, n50_contigs, contig_count, gtdb_genome_representative, checkm2_completeness, checkm2_contamination, ssu_length, ssu_query_id) %>% 
  mutate(
    closest_genome_reference = str_remove(gtdb_genome_representative, "^[A-Z]+_")
  ) %>% 
  rename(
    completeness_checkm2 = checkm2_completeness,
    contamination_checkm2 = checkm2_contamination,
    contig_n50 = n50_contigs,
    n_contigs = contig_count
  )
  

gtdb_compare = gtdb_mimag_stats_distinct %>%
  select(bin, classification:species, mimag, genome_size, contig_n50, n_contigs, closest_genome_reference, completeness_checkm2, contamination_checkm2) %>% #n_transposons
  left_join(gtdb_rep_sel, by = "closest_genome_reference") %>% 
  left_join(gtdb_redun %>% rename(classification = gtdb_taxonomy), by = "classification") %>% 
  # left_join(gtdb_me_transposons, by = "closest_genome_reference") %>% 
  # mutate(
  #   n_transposons.x = ifelse(is.na(n_transposons.x), 0, n_transposons.x),
  #   n_transposons.y = ifelse(is.na(n_transposons.y), 0, n_transposons.y)
  # ) %>% 
  filter(gtdb_type_designation_ncbi_taxa == "not type material") %>% 
  group_by(bin) %>% 
  filter(n() == 1) %>% 
  ungroup()



library(tidyr)
library(stringr)
library(scales)
library(patchwork)

df <- fread("./data/fecal_studies.csv")

palette_32 <- c(
  "#D5695D",
  "#C1534B",
  "#E38A71",
  "#B75C44",
  "#E29C74",
  "#C97A63",
  "#708238", 
  "#5D6C30",
  "#8D9B4A",
  "#A4B86B",
  "#6B8740",
  "#4C5F2A",
  "#628A8D",
  "#496F72",
  "#7FA4A7",
  "#4F6F88",
  "#6E8EAA",
  "#5A748C",
  "#D0A46B",
  "#B6895A",
  "#9E6B4F",
  "#C9B39E",
  "#9BA7A5",
  "#7C8B89",
  "#8A6A7A",
  "#A27E8E",
  "#3F5A4F",
  "#7F977F",
  "#53798C",
  "#94B0B5",
  "#C47C7D",
  "#6A7A4F"
)

p <- ggplot(df, aes(x = samples, y = total_Gbp, color = year_author, shape = technology)) +
  geom_line(aes(group = year_author)) +
  geom_point(alpha = 0.6, size = 4) +
  geom_segment(aes(x = 200, xend = 230, y = 3600, yend = 3750),
               arrow = arrow(length = unit(.15, 'cm')),
               color = "black") +
  annotate("text", y = 3600, x = 170, label = "This study",
           size = 8/.pt, fontface = "bold", family = "Times New Roman") +
  labs(
    x = "Total samples in study",
    y = "Total Gbp in study",
    color = "Study",
    shape = "Technology",
    title = "A) Fecal Long-Read Microbiome Studies",
    subtitle = NULL
    # subtitle = "Comparison across studies and sequencing technologies"
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0.4, 1, 100, 500, 1000, 3000),
    limits = c(0.4, NA),
    labels = scales::label_number()
  ) +
  guides(color = "none") +
  scale_color_manual(values = palette_32) +
  # scale_color_manual(values = scales::hue_pal()(length(unique(df$year_author)))) +
  theme(
    panel.background = element_rect(fill="grey97"),
    panel.grid.major = element_line(color = "grey85"), 
    panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
    axis.ticks.x = element_line(color = "black"), 
    axis.ticks.y = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = c(0.80, 0.10), 
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(color = "black", size = 8, family = "Times New Roman"),
    axis.text.y = element_text(color = "black", size = 8, family = "Times New Roman"), 
    axis.title.x = element_text(size = 8, family = "Times New Roman"),
    axis.title.y = element_text(size = 8, family = "Times New Roman"),
    plot.title = element_text(family = "Times New Roman", face = "bold", size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )


p

library(data.table)
library(ggplot2)

df <- fread("./data/Table_S6_LRMG_meta.csv")

palette_32 <- c(
  "#D5695D", "#C1534B", "#E38A71", "#B75C44", "#E29C74", "#C97A63",
  "#708238", "#5D6C30", "#8D9B4A", "#A4B86B", "#6B8740", "#4C5F2A",
  "#628A8D", "#496F72", "#7FA4A7", "#4F6F88", "#6E8EAA", "#5A748C",
  "#D0A46B", "#B6895A", "#9E6B4F", "#C9B39E", "#9BA7A5", "#7C8B89",
  "#8A6A7A", "#A27E8E", "#3F5A4F", "#7F977F", "#53798C", "#94B0B5",
  "#C47C7D", "#6A7A4F"
)

# Main plot with linear scale
p_main <- ggplot(df, aes(x = samples, y = total_Gbp, color = year_author, shape = technology)) +
  geom_line(aes(group = year_author)) +
  geom_point(alpha = 0.6, size = 4) +
  geom_segment(aes(x = 200, y = 3600, xend = 235, yend = 3750),
               arrow = arrow(length = unit(.15, 'cm')),
               color = "black") +
  annotate("text", y = 3600, x = 170, label = "This study",
           size = 4, fontface = "bold") +
  labs(
    x = "Total samples in study",
    y = "Total Gbp in study",
    color = "Study",
    shape = "Technology",
    title = "Fecal Long-Read Microbiome Studies"
  ) +
  guides(color = "none") +
  scale_color_manual(values = palette_32) +
  theme(
    panel.background = element_rect(fill = "grey97"),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey85", linetype = "dotted"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = "right",
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(size = 12)
  )

# Inset plot - zoomed to x=(0,50), y=(0,200)
p_inset <- ggplot(df, aes(x = samples, y = total_Gbp, color = year_author, shape = technology)) +
  geom_line(aes(group = year_author)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = palette_32) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 200)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 7),
    axis.title = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  )

# Combine plots
library(grid)
print(p_main)
print(p_inset, vp = viewport(width = 0.35, height = 0.65, x = 0.50, y = 0.4))

##################


gtdb_compare_study <- gtdb_compare %>%
  filter(mimag == "HQ") %>%
  select(bin, n_contigs.x, contig_n50.x) %>% #n_transposons.x
  pivot_longer(-bin) #%>% 
  # filter(!(str_detect(name, "n_transposons") & value == 0))

gtdb_compare_gtdb <- gtdb_compare %>%
  filter(mimag == "HQ") %>%
  select(gtdb_genome_representative, n_contigs.y, contig_n50.y) %>% #n_transposons.y
  rename(bin = gtdb_genome_representative) %>%
  distinct(bin, .keep_all = T) %>%
  pivot_longer(-bin) #%>% 
  # filter(!(str_detect(name, "n_transposons") & value == 0))

gtdb_compare_combined <- bind_rows(gtdb_compare_study, gtdb_compare_gtdb) %>%
  mutate(
    group = ifelse(str_detect(name, ".x"), "Study MAGs", "GTDB Ref."),
    metric = str_remove(name, "\\.[xy]"),
    metric = case_when(
      metric == "n_contigs" ~ "Contigs (#)",
      metric == "contig_n50" ~ "Contig N50",
      metric == "n_transposons" ~ "Mobile Elements (#)",
      TRUE ~ metric
    )
  )

# Calculate medians for annotation (before log transformation)
factor_stats <- gtdb_compare_combined %>%
  group_by(group, metric) %>%
  summarise(
    iqr = IQR(value),
    value = round(median(value)),
    .groups = "drop"
  )
total_samples <- gtdb_compare_combined %>% 
  distinct(bin, group) %>% 
  group_by(group) %>% 
  summarise(n = n()) %>% 
  mutate(metric = "Contig N50", 
         n = case_when(
           group == "GTDB Ref." ~paste0("MAGs:\n", n),
           T ~paste0("\n", n)
         ))

p_compare <- gtdb_compare_combined %>%
  ggplot(aes(x = group, y = value, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.8) +
  # geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_quasirandom(alpha = 0.4, ) +
  geom_label(
    data = factor_stats,
    aes(label = value),
    # vjust = 1.5, hjust = 0.5,
    color = "black",
    size = 7/.pt, fontface = "bold", family = "Times New Roman"
  ) +
  geom_text(
    data = total_samples,
    aes(label = n, y = 100),
    # vjust = 1.5, hjust = 0.5,
    color = "black",
    size = 7/.pt, fontface = "bold", family = "Times New Roman"
  ) +
  facet_wrap(~metric, scales = "free", ncol = 2) +
  scale_color_manual(
    values = c("Study MAGs" = "#4f6e9a", "GTDB Ref." = "#9a4f6e")
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
    limits = c(1, NA),
    labels = scales::label_number()
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Dataset",
    title = "B) High Quality MAGs vs Non-typestrain GTDB References",
  ) +
  theme(
    panel.background = element_rect(fill="grey97"),
    panel.grid.major = element_line(color = "grey85"), 
    panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
    axis.ticks.x = element_line(color = "black"), 
    axis.ticks.y = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = "none", 
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(color = "black", size = 8, family = "Times New Roman"),
    axis.text.y = element_text(color = "black", size = 8, family = "Times New Roman"), 
    axis.title.x = element_text(size = 8, family = "Times New Roman"),
    axis.title.y = element_text(size = 8, family = "Times New Roman"),
    strip.text = element_text(size = 8, family = "Times New Roman"), 
    plot.title = element_text(family = "Times New Roman", face = "bold", size = 10)
#    panel.background = element_rect(fill="grey97"),
#    panel.grid.major = element_line(color = "grey85"), 
#    panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
#    axis.ticks.x = element_line(color = "black"), 
#    axis.ticks.y = element_line(color = "black"),
#    axis.line = element_line(color = "black", linewidth = 0.1),
#    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
#    #panel.border = element_blank(),
#    strip.text = element_text(face = "bold"),
#    plot.background = element_rect(fill = "transparent", color = "transparent"),
#    legend.position = "none", 
#    axis.text.x = element_text(color = "black", size = 12),
#    axis.text.y = element_text(color = "black", size = 11), 
#    axis.title.y = element_text(size = 12)
  )


p_compare

##########################
SSU_ARC = 900
SSU_BAC = 1300

gtdb_ssu_qc <- gtdb_ssu %>% 
  mutate(
    ssu_len = str_extract(classification, "ssu_len=(\\d+)"),
    ssu_len = as.numeric(str_remove(ssu_len, "ssu_len=")),
    classification = str_remove(classification, " \\[.*")
  ) %>% 
  select(gtdb_genome_reference_id, ssu_len, classification) %>% 
  arrange(gtdb_genome_reference_id)

gtdb_ssu_qc %>% 
  separate(classification, sep = ";", into = c("domain", "phylum", "class", "order", "family", "genus", "species"), remove = F) %>% 
  group_by(domain) %>% 
  summarise(
    median = median(ssu_len),
    q25 = quantile(ssu_len, 0.25),
    q75 = quantile(ssu_len, 0.75),
    pct_above_1400 = mean(ssu_len >= 1400)
  )

gtdb_ssu_qc_agg <- gtdb_ssu_qc %>% 
  filter((str_detect(classification, "Archaea") & ssu_len >= SSU_ARC) | (str_detect(classification, "Bacteria") & ssu_len >= SSU_BAC)) %>% 
  distinct(gtdb_genome_reference_id, .keep_all = T) %>% 
  group_by(classification) %>% 
  summarise(
    n_ssu_cluster = n(),
    median_len = median(ssu_len)
  ) %>% 
  ungroup()

# 
# gtdb_mimag_stats_distinct %>% 
#   filter(species == "s__") %>% 
#   select(bin) %>% 
#   mutate(file = paste0("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/mmlong2_lite_collected_bins/",bin,".fa")) %>% 
#   relocate(file) %>% 
#   write_delim(file = "batchfile.tsv", col_names = F, delim ="\t")

gtdb_mimag_stats_distinct %>% 
  filter(closest_genome_reference != "N/A") %>% 
  group_by(domain, mimag) %>% 
  summarise(
    median = median(ssu_16s_len),
    q0 = min(ssu_16s_len),
    q25 = quantile(ssu_16s_len, 0.25),
    q75 = quantile(ssu_16s_len, 0.75),
    pct_above_1400 = mean(ssu_16s_len >= 1400)
  )


study_ssu_hqmq <- gtdb_mimag_stats_distinct %>% 
  filter(closest_genome_reference != "N/A") %>% 
  filter((str_detect(classification, "Archaea") & ssu_16s_len >= SSU_ARC) | (str_detect(classification, "Bacteria") & ssu_16s_len >= SSU_BAC)) %>% 
  group_by(classification) %>% 
  summarise(
    n_ssu_study = n(),
    median_len_study = median(ssu_16s_len)
  ) %>% 
  ungroup() %>% 
  mutate(
    group = "HQ/MQ MAGs"
  ) %>% 
  left_join(gtdb_ssu_qc_agg, by = "classification") %>% 
  left_join(gtdb_redun %>% rename(classification = gtdb_taxonomy)) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(tax_level = "Species")
  

study_ssu_hq <- gtdb_mimag_stats_distinct %>% 
  filter(closest_genome_reference != "N/A", mimag == "HQ") %>% 
  filter((str_detect(classification, "Archaea") & ssu_16s_len >= 900) | (str_detect(classification, "Bacteria") & ssu_16s_len >= 1200)) %>% 
  group_by(classification) %>% 
  summarise(
    n_ssu_study = n(),
    median_len_study = median(ssu_16s_len)
  ) %>% 
  ungroup() %>% 
  mutate(
    group = "HQ MAGs"
  ) %>% 
  left_join(gtdb_ssu_qc_agg, by = "classification") %>% 
  left_join(gtdb_redun %>% rename(classification = gtdb_taxonomy)) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(tax_level = "Species")

study_ssu <- bind_rows(study_ssu_hq, study_ssu_hqmq)


gtdb_ssu_qc_agg_genus <- gtdb_ssu_qc %>%
  filter((str_detect(classification, "Archaea") & ssu_len >= SSU_ARC) | (str_detect(classification, "Bacteria") & ssu_len >= SSU_BAC)) %>%
  mutate(classification = str_remove(classification, "s\\_\\_.*")) %>%
  # distinct(gtdb_genome_reference_id, .keep_all = T) %>%
  group_by(classification) %>%
  summarise(
    n_ssu_cluster = n(),
    median_len = median(ssu_len)
  )

gtdb_redun_genus <- gtdb_redun %>% 
  rename(classification = gtdb_taxonomy) %>% 
  mutate(classification = str_remove(classification, "s\\_\\_.*")) %>%
  group_by(classification) %>%
  summarise(
    cluster_size = sum(cluster_size)
  ) %>% 
  ungroup() %>% 
  mutate(is_singleton = ifelse(cluster_size == 1, T, F))

# HQ/MQ at genus level
study_ssu_hqmq_genus <- gtdb_mimag_stats_distinct %>%
  filter(closest_genome_reference != "N/A") %>%
  filter((str_detect(classification, "Archaea") & ssu_16s_len >= SSU_ARC) | (str_detect(classification, "Bacteria") & ssu_16s_len >= SSU_BAC)) %>%
  mutate(classification = str_remove(classification, "s\\_\\_.*")) %>%
  group_by(classification) %>%
  summarise(
    n_ssu_study = n(),
    median_len_study = median(ssu_16s_len)
  ) %>%
  ungroup() %>%
  mutate(group = "HQ/MQ MAGs") %>%
  left_join(gtdb_ssu_qc_agg_genus, by = "classification") %>%
  left_join(gtdb_redun_genus) %>% 
  mutate_all(~replace(., is.na(.), 0))%>% 
  mutate(tax_level = "Genus")

# HQ at genus level
study_ssu_hq_genus <- gtdb_mimag_stats_distinct %>%
  filter(closest_genome_reference != "N/A", mimag == "HQ") %>%
  filter((str_detect(classification, "Archaea") & ssu_16s_len >= SSU_ARC) | (str_detect(classification, "Bacteria") & ssu_16s_len >= SSU_BAC)) %>%
  mutate(classification = str_remove(classification, "s\\_\\_.*")) %>%
  group_by(classification) %>%
  summarise(
    n_ssu_study = n(),
    median_len_study = median(ssu_16s_len)
  ) %>%
  ungroup() %>%
  mutate(group = "HQ MAGs") %>%
  left_join(gtdb_ssu_qc_agg_genus, by = "classification") %>%
  left_join(gtdb_redun_genus) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(tax_level = "Genus")




library(ggalluvial)
study_ssu_hqmq_combined <- bind_rows(study_ssu_hqmq, study_ssu_hqmq_genus) %>% 
  mutate(
    before = case_when(
      n_ssu_cluster == 0 ~ "No SSU",
      n_ssu_cluster == 1 ~ "1 SSU",
      n_ssu_cluster == 2 ~ "2 SSUs",
      TRUE ~ "3+ SSUs"
    ),
    after = case_when(
      n_ssu_cluster + n_ssu_study == 0 ~ "No SSU",
      n_ssu_cluster + n_ssu_study == 1 ~ "1 SSU",
      n_ssu_cluster + n_ssu_study == 2 ~ "2 SSUs",
      TRUE ~ "3+ SSUs"
    ),
    before = factor(before, levels = c("No SSU", "1 SSU", "2 SSUs", "3+ SSUs")),
    after = factor(after, levels = c("No SSU", "1 SSU", "2 SSUs", "3+ SSUs"))
  ) %>% 
  group_by(tax_level) %>% 
  count(before, after)

study_ssu_hq_combined <- bind_rows(study_ssu_hq, study_ssu_hq_genus) %>% 
  mutate(
    before = case_when(
      n_ssu_cluster == 0 ~ "No SSU",
      n_ssu_cluster == 1 ~ "1 SSU",
      n_ssu_cluster == 2 ~ "2 SSUs",
      TRUE ~ "3+ SSUs"
    ),
    after = case_when(
      n_ssu_cluster + n_ssu_study == 0 ~ "No SSU",
      n_ssu_cluster + n_ssu_study == 1 ~ "1 SSU",
      n_ssu_cluster + n_ssu_study == 2 ~ "2 SSUs",
      TRUE ~ "3+ SSUs"
    ),
    before = factor(before, levels = c("No SSU", "1 SSU", "2 SSUs", "3+ SSUs")),
    after = factor(after, levels = c("No SSU", "1 SSU", "2 SSUs", "3+ SSUs"))
  ) %>% 
  group_by(tax_level) %>% 
  count(before, after)

# study_ssu_hq_combined_total <- study_ssu_hq_combined %>% 
#   group_by(tax_level) %>% 
#   summarise(n = sum(n)) 

p_ssu_hq <- study_ssu_hq_combined %>% 
  ggplot(aes(y = n, axis1 = before, axis2 = after)) +
  geom_alluvium(aes(fill = before), alpha = 0.7) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 1) +  # Color the strata
  # Only show counts
  geom_text(stat = "stratum", 
            aes(label = after_stat(count)),
            size = 7/.pt, fontface = "bold", color = "white", family = "Times New Roman") +
  scale_x_discrete(limits = c("GTDB R226", "GTDB+Study"), 
                   expand = c(0.15, 0.05)) +
  # geom_label(state = "stratum", data = study_ssu_hq_combined_total) +
  labs(title = "C) GTDB Clusters Gaining SSUs from HQ MAGs",
       y = NULL) +
  facet_wrap(~tax_level, scales = "free") +
  scale_fill_manual(values = c(
    "No SSU"   = "#D5695D",  # muted terracotta
    "1 SSU"    = "#D0A46B",  # warm ochre
    "2 SSUs"   = "#7FA4A7",  # soft teal-grey
    "3+ SSUs"  = "#4F6F88"   # muted steel blue
  )) +
  guides(fill = guide_legend(title = "SSU Count")) +
  theme(
    # panel.background = element_rect(fill="grey97"),
    # panel.grid.major = element_line(color = "grey85"), 
    # panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
    # axis.ticks.x = element_line(color = "black"), 
    # axis.ticks.y = element_line(color = "black"),
    # axis.line = element_line(color = "black", linewidth = 0.1),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # #panel.border = element_blank(),
    # strip.text = element_text(face = "bold"),
    # plot.background = element_rect(fill = "transparent", color = "transparent"),
    # legend.position = "bottom", 
    # axis.text.x = element_text(color = "black", size = 12),
    # axis.text.y = element_text(color = "black", size = 11), 
    # axis.title.y = element_text(size = 12)
    panel.background = element_rect(fill="grey97"),
    panel.grid.major = element_line(color = "grey85"), 
    panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
    axis.ticks.x = element_line(color = "black"), 
    axis.ticks.y = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = "none", 
    text = element_text(family = "Times New Roman"),
    strip.text = element_text(color = "black", size = 8, family = "Times New Roman"),
    axis.text.x = element_text(color = "black", size = 8, family = "Times New Roman"),
    axis.text.y = element_text(color = "black", size = 8, family = "Times New Roman"), 
    axis.title.x = element_text(size = 8, family = "Times New Roman"),
    axis.title.y = element_text(size = 8, family = "Times New Roman"),
    plot.title = element_text(size = 10, family = "Times New Roman", face = "bold")
  )

# p_ssu_hqmq <- study_ssu_hqmq_combined %>% 
#   ggplot(aes(y = n, axis1 = before, axis2 = after)) +
#   geom_alluvium(aes(fill = before), alpha = 0.7) +
#   geom_stratum(aes(fill = after_stat(stratum)), alpha = 1) +  # Color the strata
#   # Only show counts
#   geom_text(stat = "stratum", 
#             aes(label = after_stat(count)),
#             size = 3.5, fontface = "bold", color = "white") +
#   scale_x_discrete(limits = c("GTDB R226", "GTDB+Study"), 
#                    expand = c(0.15, 0.05)) +
#   # geom_label(state = "stratum", data = study_ssu_hq_combined_total) +
#   labs(title = "GTDB clusters gaining SSUs from HQ/MQ MAGs",
#        y = NULL) +
#   facet_wrap(~tax_level, scales = "free") +
#   scale_fill_manual(values = c(
#     "No SSU"   = "#D5695D",  # muted terracotta
#     "1 SSU"    = "#D0A46B",  # warm ochre
#     "2 SSUs"   = "#7FA4A7",  # soft teal-grey
#     "3+ SSUs"  = "#4F6F88"   # muted steel blue
#   )) +
#   guides(fill = guide_legend(title = "SSU Count")) +
#   theme(
#    panel.background = element_rect(fill="grey97"),
#    panel.grid.major = element_line(color = "grey85"), 
#    panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
#    axis.ticks.x = element_line(color = "black"), 
#    axis.ticks.y = element_line(color = "black"),
#    axis.line = element_line(color = "black", linewidth = 0.1),
#    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
#    #panel.border = element_blank(),
#    strip.text = element_text(face = "bold"),
#    plot.background = element_rect(fill = "transparent", color = "transparent"),
#    legend.position = "bottom", 
#    axis.text.x = element_text(color = "black", size = 12),
#    axis.text.y = element_text(color = "black", size = 11), 
#    axis.title.y = element_text(size = 12)
  #   panel.background = element_rect(fill="grey97"),
  #   panel.grid.major = element_line(color = "grey85"), 
  #   panel.grid.minor = element_line(color = "grey85", linetype = "dotted"), 
  #   axis.ticks.x = element_line(color = "black"), 
  #   axis.ticks.y = element_line(color = "black"),
  #   axis.line = element_line(color = "black", linewidth = 0.1),
  #   panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  #   #panel.border = element_blank(),
  #   plot.background = element_rect(fill = "transparent", color = "transparent"),
  #   legend.position = "none", 
  #   text = element_text(family = "Times New Roman"),
  #   axis.text.x = element_text(color = "black", size = 7, family = "Times New Roman"),
  #   axis.text.y = element_text(color = "black", size = 7, family = "Times New Roman"), 
  #   axis.title.x = element_text(size = 7, family = "Times New Roman"),
  #   axis.title.y = element_text(size = 7, family = "Times New Roman"),
  #   plot.title = element_text(family = "Times New Roman", face = "bold")
  # )

p_ssu

# combined_plot <- (p + p_compare) / (p_ssu_hqmq + p_ssu_hq) + 
#   plot_annotation(tag_levels = "A",tag_suffix = ")")
# combined_plot

combined_plot <- p + (p_compare / p_ssu_hq) + 
  plot_layout(widths = c(1, 1.5)) +  # Give right side more width
  # plot_annotation(tag_levels = "A",tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Times New Roman"))
combined_plot

ggsave(
  filename = "./figures/mag_comparison.png",  # or .pdf, .svg, etc.
  plot = combined_plot,
  width = 18.8,
  height = 15,        # adjust based on your needs
  dpi = 300,          # high resolution for publication
  units = "cm"
)




###########################
# STUDY: Study accession or unique name (alias) [study_accession]
# SAMPLE: Binned sample accession or unique name (alias) [bin]
# ASSEMBLYNAME: Unique assembly name [bin]
# ASSEMBLY_TYPE: ‘binned metagenome’ 
# COVERAGE: The estimated depth of sequencing coverage [cov]
# PROGRAM: The assembly program [program]
# PLATFORM: The sequencing platform, or comma-separated list of platforms [platform]
# MINGAPLENGTH: Minimum length of consecutive Ns to be considered a gap (optional) [mingaplength]
# MOLECULETYPE: ‘genomic DNA’, ‘genomic RNA’ or ‘viral cRNA’ (optional) [moleculetype]
# DESCRIPTION: Free text description of the genome assembly (optional) [description]
# RUN_REF: Comma separated list of run accession(s) (optional) [run_ref]

checklist <- fread("./Checklist_GSC-MIMAGS_ERC000047_1764661088249.tsv", skip = 1)
sample_accessions <- fread("./samples-2025-12-02T08_47_19.csv") %>% 
  mutate(sample = paste0("NP-", title)) %>% 
  rename(sample_id = id, sample_alias = alias) %>% 
  select(sample_id, sample)

missing_tax <- fread("./found_tax.tsv") %>% select(-taxonomy_rank_used)

gtdb_ncbi_tax <- gtdb_rep %>% 
  select(gtdb_taxonomy, ncbi_taxid, ncbi_taxonomy) %>%
  rename(classification = gtdb_taxonomy) %>% 
  mutate(
    ncbi_taxonomy = ncbi_taxonomy %>%
      str_split(";") %>%
      map_chr(~ {
        # Remove empty values and GTDB placeholders like "s__", "g__", etc.
        valid_parts <- .x[.x != "" & !str_detect(.x, "^[a-z]__$")]
        # Get the last valid part
        if (length(valid_parts) > 0) {
          tail(valid_parts, 1)
        } else {
          NA_character_
        }
      }),
    ncbi_taxonomy = str_remove(ncbi_taxonomy, "^[a-z]__")
  ) %>% 
  bind_rows(missing_tax)

sample_sheet = fread("/home/projects/cu_00014/people/sebdal/micropouch/subprojects/5.data_upload/files/hp_templates/20250930_healthypouch_sample_sheet_upload.tsv", skip = 1) %>% 
  mutate(sample = paste0("NP-", sample_title)) %>% 
  select(sample, `collection date`)


study_checklist <- gtdb_mimag_stats_distinct %>% 
  filter(mimag != "LQ") %>% 
  left_join(sample_accessions) %>% 
  left_join(sample_sheet) %>% 
  left_join(gtdb_ncbi_tax, by = "classification") %>% 
  mutate(
    sample_title = paste0("Metagenome-assembled genome from human fecal sample"),
    sample_description = paste0("This sample represents a MAG derived from the metagenomic sample ", sample_id),
    `metagenomic source` = "human gut metagenome",
    `project name` = "PRJEB98069",
    `assembly quality` = case_when(
      mimag == "HQ" ~ "Multiple fragments where gaps span repetitive regions. Presence of the 23S, 16S, and 5S rRNA genes and at least 18 tRNAs",
      mimag == "MQ" ~ "Many fragments with little to no review of assembly other than reporting of standard assembly statistics"
    ),
    `completeness software` = "CheckM2",
    `binning software` = paste0("mmlong2_lite_v",wf_v,"_",wf_mode,"flye_v2.9.4"),
    `binning parameters` = "coverage and kmer",
    `taxonomic identity marker` = "multi-marker approach",
    isolation_source = "human feces",
    `geographic location (latitude)` = 57.05,
    `geographic location (longitude)` = 9.93,
    `broad-scale environmental context` = "human-associated habitat",
    `local environmental context` = "fecal material",
    `environmental medium` = "fecal material",
    `geographic location (country and/or sea)` = "Denmark",
    `assembly software` = "flye_v2.9.4",
  ) %>% 
  rename(
    tax_id = ncbi_taxid,
    scientific_name = ncbi_taxonomy,
    sample_alias = bin,
    `sample derived from` = sample_id,
    `completeness score` = completeness_checkm2,
    `contamination score` = contamination_checkm2,
    ) %>% 
  select(
    tax_id, # taxid
    scientific_name, # scientific name
    sample_alias,
    sample_title,
    sample_description, # sample_description: "This sample represents a MAG derived from the metagenomic sample ERSXXXXX"
    `metagenomic source`,
    `sample derived from`, #ERSXXX 
    `project name`,
    `completeness score`,
    `completeness software`,
    `contamination score`,
    `binning software`,
    `assembly quality`,
    `binning parameters`,
    `taxonomic identity marker`,
    isolation_source,
    `collection date`,
    `geographic location (latitude)`,
    `geographic location (longitude)`,
    `broad-scale environmental context`,
    `local environmental context`,
    `environmental medium`,
    `geographic location (country and/or sea)`,
    `assembly software`,
    classification
  )


write_delim(study_checklist, "./files/PRJEB98069_MAG_sample_checklist.tsv", delim = "\t")

# x <- gtdb_mimag_stats_distinct %>% 
#   filter(mimag != "LQ") %>% 
#   left_join(sample_accessions) %>% 
#   left_join(sample_sheet) %>% 
#   left_join(gtdb_rep %>% select(gtdb_taxonomy, ncbi_taxid, ncbi_taxonomy) %>% rename(classification = gtdb_taxonomy), by = "classification")%>%
#   filter(is.na(ncbi_taxid)) %>% 
#   select(classification) %>% 
#   distinct() %>% 
#   write_delim("missing_tax.tsv", delim ="\t")
gtdb_mimag_stats_distinct %>% 
  filter(species == "s__") %>% 
  select(bin, classification) %>% 
  left_join(gtdb_ncbi_tax)

mimag_stats_distinct %>% 
  select(sample, wf_v, wf_mode, wf_date, bin, cov) %>% 
  rename(assemblyname=bin, coverage = cov) %>% 
  mutate(
    study = "PRJEB98069",
    assembly_type = "binned_metagenome",
    program = paste0("mmlong2_lite_v",wf_v,"_",wf_mode,"flye_v2.9.4"),
    platform = "PromethION",
    mingaplength = NA,
    moleculetype = "genomic DNA",
    description = "Fecal binned metagenome",
    run_ref = wf_date
  ) 


drep <- tibble(
  sample = list.files("/home/projects/cu_00014/people/sebdal/micropouch/subprojects/5.data_upload/output/drep")
) %>% 
  mutate(
    bin = map(.x = sample, ~(list.files(file.path("/home/projects/cu_00014/people/sebdal/micropouch/subprojects/5.data_upload/output/drep", .x, "dereplicated_genomes"))))
  ) %>% 
  unnest(bin) %>% 
  mutate(bin = str_remove(bin,  ".fa"))

mag_accessions = fread("./mag-samples-2025-12-03.csv") %>% 
  rename(bin = alias, sample = id) %>% 
  mutate(accession_num=as.numeric(str_extract(sample, "\\d+"))) %>% 
  group_by(bin) %>% 
  filter(accession_num == max(accession_num)) %>% 
  arrange(bin)

mimag_stats_distinct %>% 
  filter(mimag != "LQ") %>% 
  select(-sample) %>% 
  left_join(mag_accessions %>% select(sample, bin)) %>% 
  select(sample, wf_v, wf_mode, wf_date, bin, cov, n_contigs) %>% 
  rename(coverage = cov) %>% 
  mutate(
    study = "PRJEB98069",
    assemblyname = bin,
    assembly_type = "binned metagenome",
    program = paste0("mmlong2_lite_v", wf_v, "_", wf_mode, "_flye_v2.9.4"),
    platform = "PromethION",
    mingaplength = NA,
    moleculetype = "genomic DNA",
    description = "Fecal binned metagenome",
    run_ref = wf_date
  ) %>% 
  select(study, sample, assemblyname, assembly_type, coverage, program, 
         platform, mingaplength, moleculetype, description, run_ref, n_contigs) %>% 
  write_delim(file = "output/manifest_all.tsv", delim = "\t")



filter(mimag_stats_distinct, is_circular) %>% arrange(completeness_checkm1) %>% print(n=250)
