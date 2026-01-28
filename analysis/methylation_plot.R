library(tidyverse)

bin_quality <- read_delim("/home/projects/cu_00014/people/albmol/healthypouch/HealthyPouch_article/data/mag_qual.tsv") %>%  
  rename(sample_barcode = sample, 
         mag_quality = mimag) 

# Figure 3. Methylation plot
sample_barcode <- "PaPr00000245MP"
contig2see <- "contig_108"
sample_barcode_np <- paste0("NP-",sample_barcode) 


m <- fread(file.path("data/methylation_plot",sample_barcode,"motifs-scored-read-methylation.tsv")) %>% mutate(motif_mod = paste(motif, mod_type, mod_position, sep="_"))
bin_q <- bin_quality %>%  filter(sample_barcode == sample_barcode_np) %>%  select(bin, completeness_checkm2, contamination_checkm2, mag_quality)
c <- fread(file.path("data/methylation_plot",sample_barcode,"contig_bin.tsv"), col.names = c("contig", "bin"))
q <- fread(file.path("data/methylation_plot",sample_barcode,"bins.tsv")) %>% 
  select(bin, contig_n50, cov) %>% 
  left_join(bin_q, by="bin") %>% 
  mutate(cov = round(cov, 0))
bm <- fread(file.path("data/methylation_plot",sample_barcode,"parent-motifs.tsv"), col.names= c("motif", "mod_type", "mod_position")) %>% mutate(motif_mod = paste(motif, mod_type, mod_position, sep="_")) 
as <- fread(file.path("data/methylation_plot",sample_barcode,"assembly_info.txt")) %>% 
  dplyr::rename(contig = '#seq_name') %>%  
  select(contig, length, cov.)  %>% 
  mutate(cov = round(cov., 0))
gtdb_bac <- fread(file.path("data/methylation_plot",sample_barcode,"gtdb_bac.tsv"))

gtdb <- gtdb_bac %>%  
  filter(grepl(sample_barcode_np, user_genome)) %>%  
  filter(grepl("Escherichia coli", classification)) 

gtdb_mod <- gtdb_bac %>%  
  mutate(tax = str_extract(classification, "[^;]+$")) %>% 
  mutate(tax = str_replace(tax, "s__", "")) %>% 
  select(user_genome, tax) %>%  
  dplyr::rename(bin = user_genome)

bins_contig <- m %>%  
  left_join(c) %>% 
  mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>% 
  left_join(q) %>%
  mutate(bin = if_else(contig == contig2see, paste0(contig), bin)) %>%  
  filter(bin != "unbinned") %>%  
  distinct(bin)


meth <- m %>% 
  left_join(c) %>% 
  mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>% 
  left_join(q) %>%  
  filter(mean_read_cov >= 10) %>%
  filter(N_motif_obs >=5) %>%
  full_join(bins_contig, by="bin") %>%  
  left_join(gtdb_mod, by = "bin") %>% 
  mutate(contig_n50 = round(contig_n50/1000))


motif_weight <- meth %>%  
  filter(motif_mod %in% bm$motif_mod) %>%  
  group_by(motif_mod, bin) %>%  
  summarise(sum_motif_per_bin = sum(motif_occurences_total)) %>%  
  ungroup() %>% 
  mutate(sum_motif_per_bin = sum_motif_per_bin/1000) %>% 
  mutate(sum_motif_per_bin = if_else(
    round(sum_motif_per_bin, 0) == 0,
    round(sum_motif_per_bin, 1),
    round(sum_motif_per_bin, 0))) %>%  
  mutate(sum_motif_per_bin = paste0(sum_motif_per_bin, "k"))


contig_meth <- meth %>% 
  filter(motif_mod %in% bm$motif_mod) %>% 
  filter(contig == contig2see) %>%  
  group_by(motif_mod) %>% 
  mutate(sum_motif_per_bin = sum(motif_occurences_total)) %>%
  ungroup() %>% 
  filter(motif_mod %in% bm$motif_mod) %>%  
  mutate(group = "unbinned") %>%  
  select(-bin) %>% 
  left_join(as, by = "contig") %>%  
  dplyr::rename(bin = contig) %>%
  mutate(length = round(length/1000)) %>% 
  mutate(
    bin_clean = paste0("Plasmid<br>", bin, " (Length=", length, " kbp, Cov=", cov., ")"),
    mag_quality = "dummy",
    contig_n50 = 1234
  ) %>% 
  select(bin, motif_mod, median, group, mag_quality, contig_n50, cov., bin_clean, sum_motif_per_bin) %>% 
  mutate(sum_motif_per_bin = if_else(nchar(as.character(sum_motif_per_bin)) > 2, sum_motif_per_bin / 1000, as.numeric(sum_motif_per_bin))) %>% 
  mutate(sum_motif_per_bin = if_else(
    round(sum_motif_per_bin, 0) == 0,
    round(sum_motif_per_bin, 1),
    round(sum_motif_per_bin, 0))) %>% 
  mutate(sum_motif_per_bin = as.character(sum_motif_per_bin)) %>% 
  mutate(sum_motif_per_bin = if_else(sum_motif_per_bin < 1, paste0(sum_motif_per_bin, "k"), sum_motif_per_bin))



bin_meth <- meth %>% 
  filter(bin != "unbinned") %>%
  filter(motif_mod %in% bm$motif_mod | is.na(motif_mod)) %>% 
  group_by(bin, motif_mod) %>% 
  summarise(median = mean(median)) %>% 
  mutate(group = "binned") %>%  
  filter(!grepl("contig", bin))


## Extract the E. coli bin
ecoli_bin <- bin_meth %>%  
  filter(bin %in% gtdb$user_genome) %>% 
  left_join(
    meth %>%
      select(bin, mag_quality, contig_n50, cov, tax) %>%
      distinct(),
    by = "bin") %>%  
  mutate(bin_clean = str_replace(as.character(bin), "^.*?\\.", ""), 
         bin_clean = paste0(bin_clean = paste0("<i>", tax, "</i><br>", bin_clean, " (", mag_quality, ", N50=", contig_n50, " kbp, Cov=", cov, ")"))) 


# Bins (minus E. coli) dataframe
bin_meth_filtered <- bin_meth %>%  
  filter(!(bin %in% gtdb$user_genome))


mat <- bin_meth_filtered %>%
  select(-group) %>% 
  pivot_wider(
    names_from   = motif_mod,
    values_from  = median,
    values_fill  = list(median = 0)) %>%   
  column_to_rownames("bin") %>%
  as.matrix()


## Cluster rows and columns 
row_hc <- hclust(dist(mat))
col_hc <- hclust(dist(t(mat)))

## Set factor levels in dendrogram order
bins <- bin_meth_filtered %>%
  left_join(
    meth %>%
      select(bin, mag_quality, contig_n50, cov, tax) %>%
      distinct(),
    by = "bin"
  ) %>%  
  mutate(
    bin       = factor(bin,       levels = row_hc$labels[row_hc$order]),
    motif_mod = factor(motif_mod, levels = col_hc$labels[col_hc$order]),
    bin_clean = str_replace(as.character(bin), "^.*?\\.", ""), 
    bin_clean = paste0(bin_clean = paste0("<i>", tax, "</i>", "<br>", bin_clean, " (", mag_quality, ", N50=", contig_n50, " kbp, Cov=", cov, ")"))) %>%
  mutate(bin_clean = factor(bin_clean, levels = unique(bin_clean[order(bin)])))


## E. coli bin and plasmid dataframe
contigbin2see <- bind_rows(contig_meth, ecoli_bin) %>%
  mutate(
    motif_mod = factor(motif_mod, levels = levels(bins$motif_mod))
  ) %>%  
  mutate(order_helper = case_when(
    grepl("^contig", bin) ~ 2,  
    TRUE ~ 1
  )) %>%
  arrange(order_helper) %>%
  select(-order_helper)


new_bin_clean_levels <- c(levels(bins$bin_clean), contigbin2see$bin_clean) %>% 
  unique()

## Combine dataframes
df_final <- bind_rows(bins, contigbin2see) %>%
  mutate(bin_clean = factor(bin_clean, levels = new_bin_clean_levels)) %>% 
  left_join(motif_weight, by = c("motif_mod", "bin")) %>%
  mutate(sum_motif_per_bin.x = as.character(sum_motif_per_bin.x)) %>% 
  mutate(sum_motif_per_bin.x = if_else(is.na(sum_motif_per_bin.x), sum_motif_per_bin.y, sum_motif_per_bin.x)) %>%  
  rename(sum_motif_per_bin=sum_motif_per_bin.x) %>%  
  select(-sum_motif_per_bin.y) %>%  
  mutate(motif_mod_val = motif_mod, 
         motif_mod_val = factor(motif_mod_val, levels = levels(bins$motif_mod)))


df_plot <- expand_grid(
  motif_mod_val = unique(df_final$motif_mod_val),
  bin_clean = new_bin_clean_levels) %>% 
  left_join(df_final, by = c("motif_mod_val", "bin_clean")) %>%  
  mutate(concat = if_else(!(is.na(bin) | is.na(motif_mod)), paste0(bin, "_", motif_mod), NA_character_)) %>% 
  mutate(concat = if_else(is.na(concat), paste0(sample_barcode_np, ".", sub(" \\(.*", "", bin_clean), "_", sub(" \\(.*", "", motif_mod_val)), concat)) %>% 
  filter(motif_mod != "GACNGNNNNTTA_a_1") %>%  
  filter(!(motif_mod %in% c("CAATC_a_2", "CAAAAA_a_5"))) %>%  
  mutate(text_color = if_else(median > 0.79, "white", "black")) 


## Label modifications
bold_labels <- sapply(levels(df_final$bin_clean), function(x) {
  if (x %in% ecoli_bin$bin_clean) {
    paste0("<b>", x, "</b>")
  } else if (grepl("^Plasmid", x)) {
    paste0("<b>", x, "</b>")
  } else {
    x
  }
})

motif_labels <- c(
  "GAGCAG_a_4" = "GAGC<sub>6m</sub><b>A</b>G",
  "CTGCAG_a_4" = "CTGC<sub>6m</sub><b>A</b>G",
  "GTANNNNNNTGC_a_2" = "GT<sub>6m</sub><b>A</b>(N)<sub>6</sub>TGC",
  "RTCANNNNNNNTCRG_a_3" = "RTC<sub>6m</sub><b>A</b>(N)<sub>7</sub>TCRG",
  "GCTTGA_a_5" = "GCTTG<sub>6m</sub><b>A</b>",
  "TAANNNNNNGTCY_a_2" = "TA<sub>6m</sub><b>A</b>(N)<sub>6</sub>GTCY",
  "CTNGANG_a_4" = "CTNG<sub>6m</sub><b>A</b>NG",
  "TCCANNNNNNNTGC_a_3" = "TCC<sub>6m</sub><b>A</b>(N)<sub>7</sub>TGC",
  "GATC_a_1" = "G<sub>6m</sub><b>A</b>TC")


## Plot
fig3_plot <- df_plot  %>% 
  ggplot(aes(x = motif_mod_val, y = bin_clean, fill = median)) + 
  geom_tile(color = "gray") +
  geom_text(data = df_plot, aes(x = motif_mod_val, y = bin_clean, label = sum_motif_per_bin, color = I(text_color)), size = 3, inherit.aes = FALSE) +
  scale_fill_gradient(
    name = "Methylation level",
    low = "white",
    high = "#03396c",
    na.value = "white",
    limits = c(0, 1),              
    breaks = c(0, 0.25, 0.5, 0.75, 1.00),
    labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_x_discrete(labels = motif_labels) + 
  scale_y_discrete(limits = rev(levels(df_final$bin_clean)),
                   labels = bold_labels) +
  theme(
    axis.text.x = element_markdown(angle= 45, hjust = 1, vjust = 1, size=12), 
    axis.text.y = element_markdown(size = 12),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12, margin = margin(b = 15))) +
  coord_fixed(ratio = 1) +
  xlab("") +
  ylab("")

fig3_plot


# Supplementary figure 1. Methylation plot
sample_barcode <- "PaPr00000001MP"
contig2see <- "contig_21"
sample_barcode_np <- paste0("NP-",sample_barcode)

m <- fread(file.path("data/methylation_plot",sample_barcode,"motifs-scored-read-methylation.tsv")) %>% mutate(motif_mod = paste(motif, mod_type, mod_position, sep="_"))
bin_q <- bin_quality %>%  filter(sample_barcode == sample_barcode_np) %>%  select(bin, completeness_checkm2, contamination_checkm2, mag_quality)
c <- fread(file.path("data/methylation_plot",sample_barcode,"contig_bin.tsv"), col.names = c("contig", "bin"))
q <- fread(file.path("data/methylation_plot",sample_barcode,"bins.tsv")) %>% 
  select(bin, contig_n50, cov) %>% 
  left_join(bin_q, by="bin") %>% 
  mutate(cov = round(cov, 0))
bm <- fread(file.path("data/methylation_plot",sample_barcode,"parent-motifs.tsv"), col.names= c("motif", "mod_type", "mod_position")) %>% mutate(motif_mod = paste(motif, mod_type, mod_position, sep="_")) 
as <- fread(file.path("data/methylation_plot",sample_barcode,"assembly_info.txt")) %>% 
  dplyr::rename(contig = '#seq_name') %>%  
  select(contig, length, cov.)  %>% 
  mutate(cov = round(cov., 0))
gtdb_bac <- fread(file.path("data/methylation_plot",sample_barcode,"gtdb_bac.tsv"))


gtdb <- gtdb_bac %>%  
  filter(grepl(sample_barcode_np, user_genome)) %>%  
  filter(grepl("Escherichia coli", classification)) 

gtdb_mod <- gtdb_bac %>%  
  mutate(tax = str_extract(classification, "[^;]+$")) %>% 
  mutate(tax = str_replace(tax, "s__", "")) %>% 
  select(user_genome, tax) %>%  
  dplyr::rename(bin = user_genome)

bins_contig <- m %>%  
  left_join(c) %>% 
  mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>% 
  left_join(q) %>%
  mutate(bin = if_else(contig == contig2see, paste0(contig), bin)) %>%  
  filter(bin != "unbinned") %>%  
  distinct(bin)


meth <- m %>% 
  left_join(c) %>% 
  mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>% 
  left_join(q) %>%  
  filter(mean_read_cov >= 10) %>%
  filter(N_motif_obs >=5) %>%
  full_join(bins_contig, by="bin") %>%  
  left_join(gtdb_mod, by = "bin") %>% 
  mutate(contig_n50 = round(contig_n50/1000))


motif_weight <- meth %>%  
  filter(motif_mod %in% bm$motif_mod) %>%  
  group_by(motif_mod, bin) %>%  
  summarise(sum_motif_per_bin = sum(motif_occurences_total)) %>%  
  ungroup() %>% 
  mutate(sum_motif_per_bin = sum_motif_per_bin/1000) %>% 
  mutate(sum_motif_per_bin = if_else(
    round(sum_motif_per_bin, 0) == 0,
    round(sum_motif_per_bin, 1),
    round(sum_motif_per_bin, 0))) %>%  
  mutate(sum_motif_per_bin = paste0(sum_motif_per_bin, "k"))


contig_meth <- meth %>% 
  filter(motif_mod %in% bm$motif_mod) %>% 
  filter(contig == contig2see) %>%  
  group_by(motif_mod) %>% 
  mutate(sum_motif_per_bin = sum(motif_occurences_total)) %>%
  ungroup() %>% 
  filter(motif_mod %in% bm$motif_mod) %>%  
  mutate(group = "unbinned") %>%  
  select(-bin) %>% 
  left_join(as, by = "contig") %>%  
  dplyr::rename(bin = contig) %>%
  mutate(length = round(length/1000)) %>% 
  mutate(
    bin_clean = paste0("Plasmid<br>", bin, " (Length=", length, " kbp, Cov=", cov., ")"),
    mag_quality = "dummy",
    contig_n50 = 1234
  ) %>% 
  select(bin, motif_mod, median, group, mag_quality, contig_n50, cov., bin_clean, sum_motif_per_bin) %>%
  mutate(sum_motif_per_bin = if_else(nchar(as.character(sum_motif_per_bin)) > 2, sum_motif_per_bin / 1000, as.numeric(sum_motif_per_bin))) %>% 
  mutate(sum_motif_per_bin = if_else(
    round(sum_motif_per_bin, 0) == 0,
    round(sum_motif_per_bin, 1),
    round(sum_motif_per_bin, 0))) %>% 
  mutate(sum_motif_per_bin = as.character(sum_motif_per_bin)) %>% 
  mutate(sum_motif_per_bin = if_else(sum_motif_per_bin < 1, paste0(sum_motif_per_bin, "k"), sum_motif_per_bin))



bin_meth <- meth %>% 
  filter(bin != "unbinned") %>%
  filter(motif_mod %in% bm$motif_mod | is.na(motif_mod)) %>% 
  group_by(bin, motif_mod) %>% 
  summarise(median = mean(median)) %>% 
  mutate(group = "binned") %>%  
  filter(!grepl("contig", bin))


## Extract the E. coli bin
ecoli_bin <- bin_meth %>%  
  filter(bin %in% gtdb$user_genome) %>% 
  left_join(
    meth %>%
      select(bin, mag_quality, contig_n50, cov, tax) %>%
      distinct(),
    by = "bin") %>%  
  mutate(bin_clean = str_replace(as.character(bin), "^.*?\\.", ""), 
         bin_clean = paste0(bin_clean = paste0("<i>", tax, "</i><br>", bin_clean, " (", mag_quality, ", N50=", contig_n50, " kbp, Cov=", cov, ")"))) 


# Bins (minus E. coli) dataframe
bin_meth_filtered <- bin_meth %>%  
  filter(!(bin %in% gtdb$user_genome))


mat <- bin_meth_filtered %>%
  select(-group) %>% 
  pivot_wider(
    names_from   = motif_mod,
    values_from  = median,
    values_fill  = list(median = 0)) %>%     
  column_to_rownames("bin") %>%
  as.matrix()


## Cluster rows and columns 
row_hc <- hclust(dist(mat))
col_hc <- hclust(dist(t(mat)))

## Set factor levels in dendrogram order
bins <- bin_meth_filtered %>%
  left_join(
    meth %>%
      select(bin, mag_quality, contig_n50, cov, tax) %>%
      distinct(),
    by = "bin"
  ) %>%  
  mutate(
    bin       = factor(bin,       levels = row_hc$labels[row_hc$order]),
    motif_mod = factor(motif_mod, levels = col_hc$labels[col_hc$order]),
    bin_clean = str_replace(as.character(bin), "^.*?\\.", ""), 
    bin_clean = paste0(bin_clean = paste0("<i>", tax, "</i>", "<br>", bin_clean, " (", mag_quality, ", N50=", contig_n50, " kbp, Cov=", cov, ")"))) %>%
  mutate(bin_clean = factor(bin_clean, levels = unique(bin_clean[order(bin)])))


## E. coli bin and plasmid dataframe
contigbin2see <- bind_rows(contig_meth, ecoli_bin) %>%
  mutate(
    motif_mod = factor(motif_mod, levels = levels(bins$motif_mod))
  ) %>%  
  mutate(order_helper = case_when(
    grepl("^contig", bin) ~ 2,  
    TRUE ~ 1
  )) %>%
  arrange(order_helper) %>%
  select(-order_helper)


new_bin_clean_levels <- c(levels(bins$bin_clean), contigbin2see$bin_clean) %>% 
  unique()

## Combine dataframes
df_final <- bind_rows(bins, contigbin2see) %>%
  mutate(bin_clean = factor(bin_clean, levels = new_bin_clean_levels)) %>% 
  left_join(motif_weight, by = c("motif_mod", "bin")) %>%
  mutate(sum_motif_per_bin.x = as.character(sum_motif_per_bin.x)) %>% 
  mutate(sum_motif_per_bin.x = if_else(is.na(sum_motif_per_bin.x), sum_motif_per_bin.y, sum_motif_per_bin.x)) %>%  
  rename(sum_motif_per_bin=sum_motif_per_bin.x) %>%  
  select(-sum_motif_per_bin.y) %>%  
  mutate(motif_mod_val = motif_mod, 
         motif_mod_val = factor(motif_mod_val, levels = levels(bins$motif_mod)))


df_plot <- expand_grid(
  motif_mod_val = unique(df_final$motif_mod_val),
  bin_clean = new_bin_clean_levels) %>% 
  left_join(df_final, by = c("motif_mod_val", "bin_clean")) %>%  
  mutate(concat = if_else(!(is.na(bin) | is.na(motif_mod)), paste0(bin, "_", motif_mod), NA_character_)) %>% 
  mutate(concat = if_else(is.na(concat), paste0(sample_barcode_np, ".", sub(" \\(.*", "", bin_clean), "_", sub(" \\(.*", "", motif_mod_val)), concat)) %>% 
  filter(motif_mod != "GACNGNNNNTTA_a_1") %>%  
  filter(!(motif_mod %in% c("CAATC_a_2", "CAAAAA_a_5"))) %>%  
  mutate(text_color = if_else(median > 0.79, "white", "black")) 


## Label modifications
bold_labels <- sapply(levels(df_final$bin_clean), function(x) {
  if (x %in% ecoli_bin$bin_clean) {
    paste0("<b>", x, "</b>")
  } else if (grepl("^Plasmid", x)) {
    paste0("<b>", x, "</b>")
  } else {
    x
  }
})

motif_labels <- c(
  "GATC_m_3"       = "GAT<sub>5m</sub><b>C</b>",              
  "TCTG_m_1"       = "T<sub>5m</sub><b>C</b>TG",              
  "CCWGG_m_1"      = "C<sub>5m</sub><b>C</b>WGG",             
  "GAACNNNNNCGT_a_2" = "GA<sub>6m</sub><b>A</b>C(N)<sub>5</sub>CGT", 
  "GATC_a_1"       = "G<sub>6m</sub><b>A</b>TC"               
)

## Plot
suppl_fig1_plot <- df_plot  %>% 
  ggplot(aes(x = motif_mod_val, y = bin_clean, fill = median)) + 
  geom_tile(color = "gray") +
  geom_text(data = df_plot, aes(x = motif_mod_val, y = bin_clean, label = sum_motif_per_bin, color = I(text_color)), size = 3, inherit.aes = FALSE) +
  scale_fill_gradient(
    name = "Methylation level",
    low = "white",
    high = "#03396c",
    na.value = "white",
    limits = c(0, 1),              
    breaks = c(0, 0.25, 0.5, 0.75, 1.00),
    labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_x_discrete(labels = motif_labels) + 
  scale_y_discrete(limits = rev(levels(df_final$bin_clean)),
                   labels = bold_labels) +
  theme(
    axis.text.x = element_markdown(angle= 45, hjust = 1, vjust = 1, size=12), 
    axis.text.y = element_markdown(size = 12),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12, margin = margin(b = 15))) +
  coord_fixed(ratio = 1) +
  xlab("") +
  ylab("")

suppl_fig1_plot



# Supplementary figure 2. Methylation plot
sample_barcode <- "PaPr00000088MP"
contig2see <- "contig_3"
sample_barcode_np <- paste0("NP-",sample_barcode)

m <- fread(file.path("data/methylation_plot",sample_barcode,"motifs-scored-read-methylation.tsv")) %>% mutate(motif_mod = paste(motif, mod_type, mod_position, sep="_"))
bin_q <- bin_quality %>%  filter(sample_barcode == sample_barcode_np) %>%  select(bin, completeness_checkm2, contamination_checkm2, mag_quality)
c <- fread(file.path("data/methylation_plot",sample_barcode,"contig_bin.tsv"), col.names = c("contig", "bin"))
q <- fread(file.path("data/methylation_plot",sample_barcode,"bins.tsv")) %>% 
  select(bin, contig_n50, cov) %>% 
  left_join(bin_q, by="bin") %>% 
  mutate(cov = round(cov, 0))
bm <- fread(file.path("data/methylation_plot",sample_barcode,"parent-motifs.tsv"), col.names= c("motif", "mod_type", "mod_position")) %>% mutate(motif_mod = paste(motif, mod_type, mod_position, sep="_")) 
as <- fread(file.path("data/methylation_plot",sample_barcode,"assembly_info.txt")) %>% 
  dplyr::rename(contig = '#seq_name') %>%  
  select(contig, length, cov.)  %>% 
  mutate(cov = round(cov., 0))
gtdb_bac <- fread(file.path("data/methylation_plot",sample_barcode,"gtdb_bac.tsv"))

gtdb <- gtdb_bac %>%  
  filter(grepl(sample_barcode_np, user_genome)) %>%  
  filter(grepl("Escherichia coli", classification)) 

gtdb_mod <- gtdb_bac %>%  
  mutate(tax = str_extract(classification, "[^;]+$")) %>% 
  mutate(tax = str_replace(tax, "s__", "")) %>% 
  select(user_genome, tax) %>%  
  dplyr::rename(bin = user_genome)

bins_contig <- m %>%  
  left_join(c) %>% 
  mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>% 
  left_join(q) %>%
  mutate(bin = if_else(contig == contig2see, paste0(contig), bin)) %>%  
  filter(bin != "unbinned") %>%  
  distinct(bin)


meth <- m %>% 
  left_join(c) %>% 
  mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>% 
  left_join(q) %>%  
  filter(mean_read_cov >= 10) %>%
  filter(N_motif_obs >=5) %>%
  full_join(bins_contig, by="bin") %>%  
  left_join(gtdb_mod, by = "bin") %>% 
  mutate(contig_n50 = round(contig_n50/1000))


motif_weight <- meth %>%  
  filter(motif_mod %in% bm$motif_mod) %>%  
  group_by(motif_mod, bin) %>%  
  summarise(sum_motif_per_bin = sum(motif_occurences_total)) %>%  
  ungroup() %>% 
  mutate(sum_motif_per_bin = sum_motif_per_bin/1000) %>% 
  mutate(sum_motif_per_bin = if_else(
    round(sum_motif_per_bin, 0) == 0,
    round(sum_motif_per_bin, 1),
    round(sum_motif_per_bin, 0))) %>%  
  mutate(sum_motif_per_bin = paste0(sum_motif_per_bin, "k"))


contig_meth <- meth %>% 
  filter(motif_mod %in% bm$motif_mod) %>% 
  filter(contig == contig2see) %>%  
  group_by(motif_mod) %>% 
  mutate(sum_motif_per_bin = sum(motif_occurences_total)) %>%
  ungroup() %>% 
  filter(motif_mod %in% bm$motif_mod) %>%  
  mutate(group = "unbinned") %>%  
  select(-bin) %>% 
  left_join(as, by = "contig") %>%  
  dplyr::rename(bin = contig) %>%
  mutate(length = round(length/1000)) %>% 
  mutate(
    bin_clean = paste0("Plasmid<br>", bin, " (Length=", length, " kbp, Cov=", cov., ")"),
    mag_quality = "dummy",
    contig_n50 = 1234
  ) %>% 
  select(bin, motif_mod, median, group, mag_quality, contig_n50, cov., bin_clean, sum_motif_per_bin) %>%
  mutate(sum_motif_per_bin = if_else(nchar(as.character(sum_motif_per_bin)) > 2, sum_motif_per_bin / 1000, as.numeric(sum_motif_per_bin))) %>% 
  mutate(sum_motif_per_bin = if_else(
    round(sum_motif_per_bin, 0) == 0,
    round(sum_motif_per_bin, 1),
    round(sum_motif_per_bin, 0))) %>% 
  mutate(sum_motif_per_bin = as.character(sum_motif_per_bin)) %>% 
  mutate(sum_motif_per_bin = if_else(sum_motif_per_bin < 1, paste0(sum_motif_per_bin, "k"), sum_motif_per_bin))



bin_meth <- meth %>% 
  filter(bin != "unbinned") %>%
  filter(motif_mod %in% bm$motif_mod | is.na(motif_mod)) %>% 
  group_by(bin, motif_mod) %>% 
  summarise(median = mean(median)) %>% 
  mutate(group = "binned") %>%  
  filter(!grepl("contig", bin))


## Extract the E. coli bin
ecoli_bin <- bin_meth %>%  
  filter(bin %in% gtdb$user_genome) %>% 
  left_join(
    meth %>%
      select(bin, mag_quality, contig_n50, cov, tax) %>%
      distinct(),
    by = "bin") %>%  
  mutate(bin_clean = str_replace(as.character(bin), "^.*?\\.", ""), 
         bin_clean = paste0(bin_clean = paste0("<i>", tax, "</i><br>", bin_clean, " (", mag_quality, ", N50=", contig_n50, " kbp, Cov=", cov, ")"))) 


# Bins (minus E. coli) dataframe
bin_meth_filtered <- bin_meth %>%  
  filter(!(bin %in% gtdb$user_genome))


mat <- bin_meth_filtered %>%
  select(-group) %>% 
  pivot_wider(
    names_from   = motif_mod,
    values_from  = median,
    values_fill  = list(median = 0)) %>%     
  column_to_rownames("bin") %>%
  as.matrix()


## Cluster rows and columns 
row_hc <- hclust(dist(mat))
col_hc <- hclust(dist(t(mat)))

## Set factor levels in dendrogram order
bins <- bin_meth_filtered %>%
  left_join(
    meth %>%
      select(bin, mag_quality, contig_n50, cov, tax) %>%
      distinct(),
    by = "bin"
  ) %>%  
  mutate(
    bin       = factor(bin,       levels = row_hc$labels[row_hc$order]),
    motif_mod = factor(motif_mod, levels = col_hc$labels[col_hc$order]),
    bin_clean = str_replace(as.character(bin), "^.*?\\.", ""), 
    bin_clean = paste0(bin_clean = paste0("<i>", tax, "</i>", "<br>", bin_clean, " (", mag_quality, ", N50=", contig_n50, " kbp, Cov=", cov, ")"))) %>%
  mutate(bin_clean = factor(bin_clean, levels = unique(bin_clean[order(bin)])))


## E. coli bin and plasmid dataframe
contigbin2see <- bind_rows(contig_meth, ecoli_bin) %>%
  mutate(
    motif_mod = factor(motif_mod, levels = levels(bins$motif_mod))
  ) %>%  
  mutate(order_helper = case_when(
    grepl("^contig", bin) ~ 2,  
    TRUE ~ 1
  )) %>%
  arrange(order_helper) %>%
  select(-order_helper)


new_bin_clean_levels <- c(levels(bins$bin_clean), contigbin2see$bin_clean) %>% 
  unique()

## Combine dataframes
df_final <- bind_rows(bins, contigbin2see) %>%
  mutate(bin_clean = factor(bin_clean, levels = new_bin_clean_levels)) %>% 
  left_join(motif_weight, by = c("motif_mod", "bin")) %>%
  mutate(sum_motif_per_bin.x = as.character(sum_motif_per_bin.x)) %>% 
  mutate(sum_motif_per_bin.x = if_else(is.na(sum_motif_per_bin.x), sum_motif_per_bin.y, sum_motif_per_bin.x)) %>%  
  rename(sum_motif_per_bin=sum_motif_per_bin.x) %>%  
  select(-sum_motif_per_bin.y) %>%  
  mutate(motif_mod_val = motif_mod, 
         motif_mod_val = factor(motif_mod_val, levels = levels(bins$motif_mod)))


df_plot <- expand_grid(
  motif_mod_val = unique(df_final$motif_mod_val),
  bin_clean = new_bin_clean_levels) %>% 
  left_join(df_final, by = c("motif_mod_val", "bin_clean")) %>%  
  mutate(concat = if_else(!(is.na(bin) | is.na(motif_mod)), paste0(bin, "_", motif_mod), NA_character_)) %>% 
  mutate(concat = if_else(is.na(concat), paste0(sample_barcode_np, ".", sub(" \\(.*", "", bin_clean), "_", sub(" \\(.*", "", motif_mod_val)), concat)) %>% 
  filter(motif_mod != "GACNGNNNNTTA_a_1") %>%  
  filter(!(motif_mod %in% c("CAATC_a_2", "CAAAAA_a_5"))) %>%  
  mutate(text_color = if_else(median > 0.79, "white", "black")) 


## Label modifications
bold_labels <- sapply(levels(df_final$bin_clean), function(x) {
  if (x %in% ecoli_bin$bin_clean) {
    paste0("<b>", x, "</b>")
  } else if (grepl("^Plasmid", x)) {
    paste0("<b>", x, "</b>")
  } else {
    x
  }
})

motif_labels <- c(
  "GATC_a_1"        = "G<sub>6m</sub><b>A</b>TC",                 
  "CCWGG_m_1"       = "C<sub>5m</sub><b>C</b>WGG",                
  "GATATC_a_1"      = "G<sub>6m</sub><b>A</b>TATC",               
  "ANGTATAC_a_6"    = "ANGTAT<sub>6m</sub><b>A</b>C",
  "ATGNNNNNGCT_a_0" = "<sub>6m</sub><b>A</b>TG(N)<sub>5</sub>GCT" 
)

## Plot
suppl_fig2_plot <- df_plot  %>% 
  ggplot(aes(x = motif_mod_val, y = bin_clean, fill = median)) + 
  geom_tile(color = "gray") +
  geom_text(data = df_plot, aes(x = motif_mod_val, y = bin_clean, label = sum_motif_per_bin, color = I(text_color)), size = 3, inherit.aes = FALSE) +
  scale_fill_gradient(
    name = "Methylation level",
    low = "white",
    high = "#03396c",
    na.value = "white",
    limits = c(0, 1),              
    breaks = c(0, 0.25, 0.5, 0.75, 1.00),
    labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_x_discrete(labels = motif_labels) + 
  scale_y_discrete(limits = rev(levels(df_final$bin_clean)),
                   labels = bold_labels) +
  theme(
    axis.text.x = element_markdown(angle= 45, hjust = 1, vjust = 1, size=12), 
    axis.text.y = element_markdown(size = 12),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12, margin = margin(b = 15))) +
  coord_fixed(ratio = 1) +
  xlab("") +
  ylab("")

suppl_fig2_plot


