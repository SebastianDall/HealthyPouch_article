library(tidyverse)
library(data.table)
# library(jsonlite)
# library(mplibrary)
# library(ggrepel)
# library(vegan)
# library(ggpubr)



# metadata <- read_delim("data/metadata/metadata.csv")
# included <- read_delim("data/metadata/included_participants_hp.csv")
# fmtp <- fread("data/metadata/redcap/FMTcohortpatient_DATA_2025-02-20_1250.csv") %>% 
#   rename(id = patient_id, sex = gender) %>% 
#   mutate(
#     id = as.character(id),
#     id = case_when(
#       nchar(id) == 1 ~paste0("pt00", id),
#       nchar(id) == 2 ~paste0("pt0", id),
#     )
#   ) %>% 
#   filter(!is.na(sex)) %>% 
#   left_join(
#     included %>% filter(str_detect(sample_barcode, "FMTP"))
#   ) %>% 
#   filter(!is.na(sample_barcode))
#
# don <- fread("/home/projects/cu_00014/people/sebdal/micropouch/subprojects/2.metadata/data/metadata/redcap/FMTdonor_DATA_2025-06-18_1149.csv") %>% 
#   rename(id = donor_id, sex = gender) %>% 
#   mutate(
#     id = as.character(id),
#     id = case_when(
#       nchar(id) == 1 ~paste0("donor00", id),
#       nchar(id) == 2 ~paste0("donor0", id),
#     )
#   ) %>% 
#   filter(!is.na(sex)) %>% 
#   left_join(
#     included %>% filter(str_detect(sample_barcode, "DON"))
#   ) %>% 
#   filter(!is.na(sample_barcode)) %>% 
#   relocate(sample_barcode, .after  = id)
#
#
assign_group <- function(df) {
  df2 <- df %>% 
    mutate(
      study_group = case_when(
        str_detect(str_to_lower(sample_barcode), "rask") ~ "normal_gut",
        str_detect(str_to_lower(sample_barcode), "do") ~ "normal_gut",

        str_detect(str_to_lower(sample_barcode), "asym") ~ "healthy_pouch",

        str_detect(str_to_lower(sample_barcode), "pa") ~ "sick_pouch",
        str_detect(str_to_lower(sample_barcode), "pouc") ~ "sick_pouch",
        str_detect(str_to_lower(sample_barcode), "fmt") ~ "sick_pouch",

        str_detect(str_to_lower(sample_barcode), "uc") ~ "uc",
      )
    )
  stopifnot(!all(is.na(df2$study_group)))
  return(df2)
}
#
# included %>% 
#   assign_group() %>% 
#   group_by(study_group) %>% 
#   summarise(n())
#
#
# ########## AGE GENDER
# do <- metadata %>%
#   filter(str_detect(id, "do")) %>%
#   select(id, sex, age) %>% 
#   filter(!is.na(sex)) %>% 
#   distinct()
#
# do_inc <- included %>% filter(str_detect(sample_barcode, "DoSc")) %>% 
#   left_join(do)
#
# age_gender <- included %>% 
#   left_join(metadata %>% select(id, sample_barcode, sex, age)) %>% 
#   distinct(sample_barcode, .keep_all = T) %>% 
#   filter(!str_detect(sample_barcode, "FMTP")) %>% 
#   filter(!str_detect(sample_barcode, "DoSc")) %>%
#   filter(!str_detect(sample_barcode, "DON")) %>%
#   bind_rows(fmtp %>% select(id, sample_barcode, age, sex)) %>% 
#   bind_rows(don %>% select(id, sample_barcode, age, sex)) %>% 
#   bind_rows(do_inc) %>% 
#   assign_group()
#
#
#
# age_gender %>% 
#   # filter(!is.na(sex)) %>% 
#   filter(age < 100, age > 1) %>% 
#   group_by(study_group) %>% 
#   summarise(
#     m = mean(age),
#     min = min(age),
#     max = max(age)
#   )
#
# age_gender %>% 
#   filter(!is.na(sex)) %>% 
#   mutate(sex = ifelse(sex == 2, 0, sex)) %>% 
#   group_by(study_group, sex) %>% 
#   summarise(
#     n()
#   )
#
# metadata_stool_fre <- metadata %>%
#   filter(sample_barcode %in% included$sample_barcode) %>% 
#   distinct(sample_barcode, .keep_all = T) %>% 
#   select(id, sample_barcode, stool_fre, bristol, contains("pdai")) %>% 
#   assign_group() %>% 
#   filter(!is.na(stool_fre)) %>% 
#   mutate(stool_fre = str_remove(stool_fre, "\\+")) %>% 
#   mutate(stool_fre_max = if_else(str_detect(stool_fre, "-"), 
#                                  map_chr(str_split(stool_fre, "-"), ~ as.character(max(as.numeric(.)))), 
#                                  as.character(stool_fre))) %>%  
#   mutate(stool_fre_min = if_else(str_detect(stool_fre, "-"),
#                                  map_chr(str_split(stool_fre, "-"), ~ as.character(min(as.numeric(.)))), 
#                                  as.character(stool_fre))) %>% 
#   mutate(stool_fre_min = as.numeric(stool_fre_min)) %>% 
#   mutate(stool_fre_max = as.numeric(stool_fre_max)) %>% 
#   mutate(stool_fre_mod = rowMeans(select(., stool_fre_min, stool_fre_max), na.rm = TRUE)) %>% 
#   mutate(
#     stool_fre_category = case_when(
#       stool_fre_mod >=   0 & stool_fre_mod <=  3 ~ "≤3",
#       stool_fre_mod >    3 & stool_fre_mod <=  5 ~ "≤5",
#       stool_fre_mod >    5 & stool_fre_mod <=  8 ~ "≤8",
#       stool_fre_mod >    8 & stool_fre_mod <= 11 ~ "≤11",
#       stool_fre_mod >   11                       ~ "≥12",
#       TRUE                                       ~ NA_character_
#     ),
#     stool_fre_category = factor(
#       stool_fre_category,
#       levels = c("≤3", "≤5", "≤8", "≤11", "≥12"),
#       ordered = TRUE
#     )
#   )
#
# metadata_stool_fre %>% 
#   group_by(study_group) %>% 
#   summarise(
#     m = mean(stool_fre_mod),
#     mi = min(stool_fre_min),
#     ma = max(stool_fre_max),
#     med = median(stool_fre_mod)
#   )
#
# metadata_stool_fre %>% 
#   mutate(study_group = factor(study_group, levels = c("sick_pouch", "healthy_pouch", "uc", "normal_gut"))) %>% 
#   ggplot(aes(x = study_group, y = stool_fre_mod)) + 
#   geom_boxplot(outliers = F) +
#   geom_point(alpha = 0.6, position = position_jitter())
#
# my_comparisons <- list(
#   c("sick_pouch",   "healthy_pouch"),
#   c("sick_pouch",   "uc"),
#   c("sick_pouch",   "normal_gut"),
#   c("healthy_pouch","uc"),
#   c("healthy_pouch","normal_gut"),
#   c("uc",           "normal_gut")
# )
#
# metadata_stool_fre %>% 
#   mutate(study_group = factor(study_group, 
#                               levels = c("sick_pouch", "healthy_pouch", "uc", "normal_gut"))) %>% 
#   ggplot(aes(x = study_group, y = stool_fre_mod)) + 
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(alpha = 0.6, width = 0.2) +
#   stat_compare_means(
#     comparisons = my_comparisons,
#     method      = "t.test",
#     label       = "p.signif",    # or "p.format" for exact p-values
#     symnum.args = list(
#       cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#       symbols   = c("p < 0.001", "p < 0.01", "p < 0.05", "ns")
#     ),
#     tip.length  = 0.02,
#     size        = 3
#   ) +
#   theme_minimal() +
#   labs(x = NULL,
#        y = "Stool Frequency",
#        title = "Pairwise t-tests between Study Groups")
#
#
# x <- metadata_stool_fre %>% 
#   filter(study_group == "normal_gut") %>% 
#   arrange(desc(stool_fre_max))


######### ILLL ########
# ill <- included %>% 
#   mutate(
#     read_count_files = map(sample_barcode, ~list.files(file.path("data/metaphlan4/", .x, "/"), pattern = "ReadCount", full.names = T))
#   ) %>% 
#   unnest(read_count_files) %>% 
#   mutate(
#     ill_read_count = map(
#       read_count_files, ~(
#         fread(
#           .x,
#           header = F,
#           col.names = c("sample_barcode2", "lib", "num_reads")
#         )
#       )
#     )
#   ) %>% 
#   unnest(ill_read_count) %>% 
#   filter(sample_barcode == sample_barcode2) %>% 
#   select(id, sample_barcode, num_reads) %>% 
#   assign_group()
# 
# 
# 
# ill_depth <- ill %>% 
#   group_by(sample_barcode) %>% 
#   summarise(num_reads = sum(num_reads)) %>% 
#   ungroup() %>% 
#   mutate(
#     bp = num_reads * 151
#   ) %>% 
#   # filter(num_reads >= 2500000) %>% 
#   assign_group()

# write_delim(ill_depth, file = "data/illumina_data.tsv", delim = "\t")

ill_depth <- fread("data/illumina_data.tsv") %>%
    filter(sample_barcode != "FMTP00000004MP")

ill_depth %>% 
  distinct(sample_barcode, .keep_all = T) %>% 
  group_by(study_group) %>% 
  summarise(n())


ill_depth %>% 
  pivot_longer(cols = -c(sample_barcode, study_group)) %>% 
  group_by(study_group, name) %>% 
  summarise(
    mean = mean(value),
    median = median(value)
  )%>% 
  arrange(name) %>% 
  mutate(
    mean = case_when(
      name == "bp" ~ round(mean / 10^9, 1),
      name == "num_reads" ~ round(mean / 10^6, 1),
    )
  )
ill_depth %>% 
  summarise(
    sum(bp / 10^9)
  )
############## NP ####################

# load_json <- function(path) {
#   # if the file doesn’t exist, return an empty list
#   if (!file.exists(path)) {
#     return(list())
#   }
#   # otherwise try to parse it, but on any error also return empty list
#   tryCatch(
#     fromJSON(path),
#     error = function(e) {
#       warning(sprintf("Could not parse JSON at %s: %s", path, e$message))
#       list()
#     }
#   )
# }


# load file if exists
# load_file <- function(file) {
#   if (file.exists(file)) {
#     fread(file)
#   } else {
#     tibble()
#   }
# }


# all_samples <- tibble(file = list.dirs("data/healthypouch/np_data/sample_barcodes", recursive = F)) %>% 
#   filter(!str_detect(file, "MP018")) %>% 
#   mutate(
#     sample_barcode = basename(file)
#   ) %>% 
#   mutate(
#     nanoq_files = map(sample_barcode,
#                       ~(load_json(file.path("data/healthypouch/np_data/sample_barcodes", .x, paste0(.x,".json"))))),
#   )
# 
# 
# np_data <- all_samples %>% 
#   mutate(
#     nanoq_files = map(sample_barcode,
#                       ~(load_json(file.path("data/healthypouch/np_data/sample_barcodes", .x, paste0(.x,".json"))))),
#     assembly_stats = map(
#       sample_barcode,
#       ~(load_file(file.path("output/assembly_stats/", paste0(.x, "_stats.csv"))))
#     )
#   )
# 
# d_np <- all_samples %>% 
#   filter(map_dbl(nanoq_files, ~(length(.x))) > 0) %>%
#   mutate(
#     read_n50 = map_dbl(nanoq_files, ~(.x$n50)),
#     bp = map_dbl(nanoq_files, ~(.x$bases)) / 10^9
#   ) %>% 
#   select(-c(nanoq_files))
# 
# d_np %>% 
#   summarise(
#     sumgb = sum(bp),
#     mean_gb = mean(bp),
#     n = n()
#   )
# 
# 
# # assembly_stats <- fread("output/assembly_stats/stats_combined.csv")
# 
# d_np <- np_data %>% 
#   select(-file) %>% 
#   # filter(map_dbl(nanoq_files, ~(length(.x))) > 0, map_dbl(assembly_stats, ~(length(.x))) > 0) %>%
#   mutate(
#     read_n50 = map_dbl(nanoq_files, ~(.x$n50)),
#     bp = map_dbl(nanoq_files, ~(.x$bases)) / 10^9
#   ) %>% 
#   assign_group() %>% 
#   select(-c(nanoq_files)) %>%
#   unnest(assembly_stats, keep_empty = TRUE) %>% 
#   mutate(
#     perc_binned = total_length_binned / total_length,
#     total_length = round(total_length / 10^9, 3)
#   )


# write_delim(d_np, "./data/nanopore_sequencing_stats.tsv", delim="\t")
d_np <- fread("./data/nanopore_sequencing_stats.tsv") %>%
    filter(sample_barcode != "FMTP00000004MP") 

d_np %>% 
  filter(!is.na(num_contigs)) %>%
  summarise(
  sumgb = sum(bp),
  mean_gb = mean(bp),
  n = n()
)

d_np %>% 
  filter(!is.na(num_contigs)) %>%
  group_by(study_group) %>% 
  summarise(n())

d_np %>% 
  pivot_longer(cols = -c(sample_barcode, study_group)) %>% 
  ggplot(aes(x = study_group, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(alpha = 0.5) + 
  facet_wrap(~name, scales = "free")
  

stats <- d_np %>% 
  filter(!is.na(num_contigs)) %>%
  pivot_longer(cols = -c(sample_barcode, study_group)) %>% 
  group_by(study_group, name) %>% 
  summarise(
    n = n(),
    mean = mean(value),
    std = sd(value),
    median = median(value),
    iqr = IQR(value),
  ) %>% 
  select(-mean, -std) %>% 
  pivot_longer(cols = c(n, median, iqr), names_to = "metric", values_to = "value") %>%
  # pivot_longer(cols = c(mean, std), names_to = "metric", values_to = "value") %>%
  pivot_wider(id_cols = study_group, names_from = c(name, metric), values_from = value)
  # print(n=30)
  # filter(study_group == "sick_pouch")


np_data_5khz <- fread("data/np_5khz_samples.txt", header = F, col.names = c("sample_barcode")) %>% 
    filter(sample_barcode != "FMTP00000004MP") %>%
    assign_group()


np_data_5khz %>% 
  group_by(study_group) %>% 
  summarise(n())




##########################################################


metaphlan <- fread("data/metaphlan4/MetaPhlAn_4.1.0_Combined_NonHuman_Subsampled_2500000_profile.txt") %>% 
  filter_taxonomy("species") %>%
  remove_NonHuman_from_colnames() %>% 
  select(clade_name, any_of(included$sample_barcode))

metaphlan_sub <- metaphlan %>% 
  column_to_rownames("clade_name")

varying_species <- apply(metaphlan_sub, 1, var) > 0
metaphlan_sub2   <- metaphlan_sub[varying_species, ]


metaphlan_hel <- decostand(t(metaphlan_sub2), method = "hellinger")


# 4. Run PCA on the transposed matrix (samples × species)
pca_res <- prcomp(metaphlan_hel, center = T, scale =F)

# 5. Build a plotting data frame
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("sample_barcode") %>% 
  left_join(metadata_stool_fre) %>% 
  assign_group() %>% 
  relocate(stool_fre:stool_fre_category, .after = sample_barcode)

percent_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

# 6. Plot PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2, shape = study_group, color = stool_fre_category)) +
  geom_point(size = 3) +
  geom_point(data = pca_df %>% filter(sample_barcode %in% x$sample_barcode), color = "red", size = 3) +
  geom_point(data = pca_df %>% filter(sample_barcode == "Asym00000007MP"), color = "orange", size = 3) +
  # geom_text_repel(aes(label = sample), max.overlaps = 20) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of MetaPhlAn Species Abundances",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  )


pca_df %>% filter(study_group == "healthy_pouch") %>% 
  mutate(pc_group = ifelse(PC2 > 0.2, "healthyhealthy", "sickhealthy")) %>% 
  ggplot(aes(x = pc_group, y = cpdai_score)) +
  # geom_boxplot() +
  geom_point(position = position_jitter(), alpha = 0.6)


metaphlan_genus <- fread("data/metaphlan4/MetaPhlAn_4.1.0_Combined_NonHuman_Subsampled_2500000_profile.txt") %>% 
  filter_taxonomy("genus") %>%
  remove_NonHuman_from_colnames() %>% 
  select(clade_name, any_of(included$sample_barcode))

m <- metadata %>% 
  filter(project == "NP" | str_detect(sample_barcode, "Asym"))

x <- metaphlan_genus %>% 
  select(clade_name, contains("Asym", ignore.case = F)) %>%
  arrange(desc(Asym00000007MP)) %>% 
  filter(clade_name == "Segatella") %>% 
  pivot_longer(-clade_name) %>% 
  pivot_wider(names_from = clade_name, values_from = value) %>% 
  summarise(mean(Segatella))

mphlan <- metaphlan_genus %>%
  pivot_longer(-clade_name, names_to = "sample_barcode") %>% 
  filter(sample_barcode %in% included$sample_barcode) %>% 
  pivot_wider(names_from = clade_name, values_from = value) %>% 
  assign_group() %>% 
  left_join(pca_df %>% select(sample_barcode, PC2, study_group)) %>% 
  left_join(metadata_stool_fre) %>% 
  relocate(study_group:stool_fre_category, .after = sample_barcode)

mphlan %>% 
  filter(study_group == "healthy_pouch") %>% 
  mutate(pc_group = ifelse(Segatella > 1, "healthyhealthy", "sickhealthy")) %>% 
  ggplot(aes(x = sample_barcode, y = Segatella, color = pc_group)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

x <- mphlan %>% 
  filter(study_group == "healthy_pouch") %>% 
  mutate(pc_group = ifelse(Segatella > 1, "healthyhealthy", "sickhealthy")) %>% 
  filter(pc_group == "healthyhealthy")


mphlan %>% 
  filter(study_group == "healthy_pouch") %>% 
  mutate(pc_group = ifelse(Segatella > 1, "Prevotella High", "Prevotella Low")) %>% 
  select(study_group, pc_group, Segatella, stool_fre_mod, cpdai_score) %>% 
  pivot_longer(-c(study_group, pc_group)) %>% 
  ggplot(aes(x = name, y = value, color = pc_group)) +
  geom_boxplot() +
  geom_point(aes(group = pc_group), position = position_jitterdodge()) +
  facet_wrap(~name, scales = "free")


sym_args <- list(
  cutpoints = c(0,    0.001, 0.01, 0.05, 1),
  symbols   = c("p<0.001", "p<0.01", "p<0.05", "ns")
)

# one pairwise comparison
my_comparison <- list(c("Prevotella Low", "Prevotella High"))

mphlan %>% 
  filter(study_group == "healthy_pouch") %>% 
  mutate(
    pc_group = ifelse(Segatella > 5, "Prevotella High", "Prevotella Low")
  ) %>% 
  select(study_group, pc_group, Segatella, stool_fre_mod, cpdai_score) %>% 
  pivot_longer(
    cols = -c(study_group, pc_group),
    names_to  = "name",
    values_to = "value"
  ) %>% 
  mutate(
    name = factor(name, levels = c("Segatella", "cpdai_score", "stool_fre_mod"), labels = c("Prevotella", "cPDAI", "Stool Freq."))
  ) %>% 
  ggplot(aes(x = pc_group, y = value, color = pc_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  stat_compare_means(
    comparisons = my_comparison,
    method      = "t.test",
    label       = "p.signif",
    symnum.args = sym_args,
    tip.length  = 0.02,
    size        = 3
  ) +
  facet_wrap(~name, scales = "free") +
  theme_minimal() +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Healthy pouch: Prevotella High vs Low"
  ) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold"),
    axis.text.x = element_blank()
  )

