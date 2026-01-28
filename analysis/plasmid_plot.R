library(tidyverse)
library(circlize)

# Figure 3. Plasmid plot
sample_barcode <- "PaPr00000245MP"

mog_results <- fread(file.path("data/plasmid_plot",sample_barcode,"mog.tsv"))
contig_arg <- fread(file.path("data/plasmid_plot",sample_barcode,"contig_arg.tsv"))

bp_to_degrees <- function(bp, plasmid_length) {
  (bp / plasmid_length) * 360 * -1}

mog <- mog_results %>%  
  select(-annotation) %>%  
  mutate(coverage = (alignment_length / refseq_length)*100) %>%  
  group_by(contig_id) %>%  
  filter(percent_identity == max(percent_identity)) %>% 
  filter(coverage == max(coverage)) %>% 
  ungroup() %>%  
  distinct(contig_id, .keep_all = T) %>% 
  dplyr::rename(function_category = `function`) %>%  
  mutate(gene_start = as.numeric(gene_start), 
         gene_end = as.numeric(gene_end)) %>%  
  mutate(gene_name = if_else(gene_name == "NA:Keyword", "NA", gene_name)) %>%  
  mutate(plasmid = "IncI1-I(Alpha)") %>%  
  full_join(contig_arg) %>%  
  filter(function_category != "phage") %>%  
  dplyr::rename(strand_og = strand) %>%  
  mutate(starnd = if_else(strand_og == "1", "+", "-")) %>%  
  mutate(function_category = str_replace(function_category, "^([a-zA-Z])", ~str_to_upper(.x))) %>% 
  relocate(gene_start:gene_end, .after = "gene_name") 


colors <- c(
  "#B39A88",  
  "#6E6436",
  "#9C9EB5",
  "#5F7F6E",
  "#2A385B" 
)

functions <- unique(mog$function_category)
palette <- setNames(colors, functions)

bp_per_char_threshold <- 140

mog <- mog %>%
  mutate(
    color = palette[function_category],
    width_bp = gene_end - gene_start,
    label_fits_inside = width_bp > nchar(gene_name) * bp_per_char_threshold,
    label_y = ifelse(label_fits_inside, 0.5, -2),
    label_facing = ifelse(label_fits_inside, "inside", "reverse.clockwise"),
    sector = "plasmid"
  ) %>% 
  mutate(text_col = if_else(toupper(color) == "#2A385B", "white", "black"))


# Add degree columns t
plasmid_length <- 90738
mog <- mog %>%
  mutate(
    deg_start = bp_to_degrees(gene_start, plasmid_length),
    deg_end = bp_to_degrees(gene_end, plasmid_length)
  )

mog_traO <- mog %>% 
  filter(gene_name == "traO") %>% 
  group_by(gene_name) %>% 
  mutate(deg_start = min(deg_start) + 360,
         deg_end = max(deg_end)) %>% 
  ungroup() %>% 
  filter(row_number() == 1)

mog <- mog %>% 
  filter(gene_name != "traO") %>% 
  bind_rows(mog_traO)


## Plot
par(mar = c(1, 1, 1, 1), family="Times New Roman")
plot(c(-1, 1), c(-1, 1), type = "n", axes = FALSE, ann = FALSE, asp = 1)

rou1=0.875
rou2 = 0.680
pad  <- 0.02

gene_rou1 <- rou1 - pad
gene_rou2 <- rou2 + pad

draw.sector(0, 360, rou1, rou2, border = "black")

## Draw each gene
for (i in 1:nrow(mog)) {
  draw.sector(
    mog$deg_start[i],
    mog$deg_end[i],
    rou1 = gene_rou1,
    rou2 = gene_rou2,
    col = mog$color[i],
    border = "black",
    
  )
}


r_inside  <- (gene_rou1 + gene_rou2) / 2 
r_outside <- rou2 - 0.06 

for (i in 1:nrow(mog)) {
  
  label <- mog$gene_name[i]
  deg <- mog$deg_end[i] - mog$deg_start[i]
  print(deg)
  
  if (str_detect(label, paste(c("tra", "rep", "mob"), collapse = "|")) & deg < -4 | label == "SHV-12") {
    
    mid_deg <- (mog$deg_start[i] + mog$deg_end[i]) / 2
    mid_rad <- mid_deg * pi / 180
    
    if (label == "SHV-12") {
      r <- r_outside - 0.1
      color <- "black"
      cex <- 1
    } else {
      r <- r_inside
      color = "black"
      cex <- 0.7
    }
    
    x <- r * cos(mid_rad)
    y <- r * sin(mid_rad)
    
    rotation <- mid_deg
    if (mid_deg < -90 && mid_deg > -270) rotation <- rotation + 180
    
    text(
      x, y, label,
      cex = cex,
      srt = rotation,
      font = 3,
      col = color,
      family = "Times New Roman",
    )
  }
}

op <- par(family = "Times New Roman", ps = 12)
on.exit(par(op), add = TRUE)
legend(x = -0.5, y=0.45, legend = names(palette), fill = palette, border = "black", cex = 0.8, bg="white", bty = "n")

dev.off()


# Supplementary figure 1. Methylation plot
sample_barcode <- "PaPr00000001MP"

mog_results <- fread(file.path("data/plasmid_plot",sample_barcode,"mog.tsv"))
contig_arg <- fread(file.path("data/plasmid_plot",sample_barcode,"contig_arg.tsv"))

mog <- mog_results %>%  
  select(-annotation) %>%  
  mutate(coverage = (alignment_length / refseq_length)*100) %>%  
  group_by(contig_id) %>%  
  filter(percent_identity == max(percent_identity)) %>% 
  filter(coverage == max(coverage)) %>% 
  ungroup() %>%  
  distinct(contig_id, .keep_all = T) %>% 
  dplyr::rename(function_category = `function`) %>%  
  mutate(gene_start = as.numeric(gene_start), 
         gene_end = as.numeric(gene_end)) %>%  
  mutate(gene_name = if_else(gene_name == "NA:Keyword", "NA", gene_name)) %>%  
  mutate(plasmid = "plasmid") %>%  
  full_join(contig_arg) %>%  
  filter(function_category != "phage") %>%  
  dplyr::rename(strand_og = strand) %>%  
  mutate(starnd = if_else(strand_og == "1", "+", "-")) %>%  
  mutate(function_category = str_replace(function_category, "^([a-zA-Z])", ~str_to_upper(.x))) %>% 
  mutate(gene_name = if_else(gene_name == "TF1-34_00113", "NA", gene_name)) %>%  
  relocate(gene_start:gene_end, .after = "gene_name")



functions <- unique(mog$function_category)
palette <- setNames(colors, functions)

bp_per_char_threshold <- 350

mog <- mog %>%
  mutate(
    color = palette[function_category],
    width_bp = gene_end - gene_start,
    label_fits_inside = width_bp > nchar(gene_name) * bp_per_char_threshold,
    label_y = if_else(width_bp > nchar(gene_name) * bp_per_char_threshold, 0.5, NA_real_),
    label_facing = if_else(width_bp > nchar(gene_name) * bp_per_char_threshold, "inside", NA_character_)
  ) %>%  
  mutate(row_id = row_number())

not_inside <- mog %>%
  filter(!label_fits_inside) %>%
  mutate(
    alt_index = row_number(),
    label_y_alt = if_else(alt_index %% 2 == 1, 2, -2),
    label_facing_alt = if_else(alt_index %% 2 == 1, "clockwise", "reverse.clockwise")
  ) %>%
  select(row_id, label_y_alt, label_facing_alt)


mog <- mog %>%
  left_join(not_inside, by = "row_id") %>%
  mutate(
    label_y = coalesce(label_y_alt, label_y),
    label_facing = coalesce(label_facing_alt, label_facing)
  ) %>%
  select(-row_id, -label_y_alt, -label_facing_alt) %>%
  mutate(sector = "plasmid")

circos.clear()

circos.par("gap.degree" = 0, "start.degree" = 1)
circos.initialize(factors = mog$sector, xlim = c(0, 119538))


circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.1, bg.border = "black",
  panel.fun = function(x, y) {
    for (i in 1:nrow(mog)) {
      start <- mog$gene_start[i]
      end <- mog$gene_end[i]
      label <- mog$gene_name[i]
      color <- mog$color[i]
      label_y <- mog$label_y[i]
      label_facing <- mog$label_facing[i]
      
      circos.rect(start, 0, end, 1, col = color, border = "black")

      if (!is.na(label) && label != "" && label != "NA") {
        circos.text((start + end) / 2, label_y, labels = label,
                    facing = label_facing, niceFacing = TRUE, cex = 1.2, font = 3)
      }
    }
  }
)

legend("center", legend = names(palette), fill = palette, border = "black", cex = 1.2)

dev.off()


# Supplementary figure 2. Methylation plot
sample_barcode <- "PaPr00000088MP"

mog_results <- fread(file.path("data/plasmid_plot",sample_barcode,"mog.tsv"))
contig_arg <- fread(file.path("data/plasmid_plot",sample_barcode,"contig_arg.tsv"))

mog <- mog_results %>%  
  select(-annotation) %>%  
  mutate(coverage = (alignment_length / refseq_length)*100) %>%  
  group_by(contig_id) %>%  
  filter(percent_identity == max(percent_identity)) %>% 
  filter(coverage == max(coverage)) %>% 
  ungroup() %>%  
  distinct(contig_id, .keep_all = T) %>% 
  dplyr::rename(function_category = `function`) %>%  
  mutate(gene_start = as.numeric(gene_start), 
         gene_end = as.numeric(gene_end)) %>%  
  mutate(gene_name = if_else(gene_name == "NA:Keyword", "NA", gene_name)) %>%  
  mutate(plasmid = "plasmid") %>%  
  full_join(contig_arg) %>%  
  filter(function_category != "phage") %>%  
  dplyr::rename(strand_og = strand) %>%  
  mutate(starnd = if_else(strand_og == "1", "+", "-")) %>%  
  mutate(function_category = str_replace(function_category, "^([a-zA-Z])", ~str_to_upper(.x))) %>% 
  mutate(gene_name = if_else(gene_name == "TF1-34_00113", "NA", gene_name)) %>%  
  relocate(gene_start:gene_end, .after = "gene_name")



functions <- unique(mog$function_category)
palette <- setNames(colors, functions)

bp_per_char_threshold <- 350

mog <- mog %>%
  mutate(
    color = palette[function_category],
    width_bp = gene_end - gene_start,
    label_fits_inside = width_bp > nchar(gene_name) * bp_per_char_threshold,
    label_y = if_else(width_bp > nchar(gene_name) * bp_per_char_threshold, 0.5, NA_real_),
    label_facing = if_else(width_bp > nchar(gene_name) * bp_per_char_threshold, "inside", NA_character_)
  ) %>%  
  mutate(row_id = row_number())

not_inside <- mog %>%
  filter(!label_fits_inside) %>%
  mutate(
    alt_index = row_number(),
    label_y_alt = if_else(alt_index %% 2 == 1, 2, -2),
    label_facing_alt = if_else(alt_index %% 2 == 1, "clockwise", "reverse.clockwise")
  ) %>%
  select(row_id, label_y_alt, label_facing_alt)


mog <- mog %>%
  left_join(not_inside, by = "row_id") %>%
  mutate(
    label_y = coalesce(label_y_alt, label_y),
    label_facing = coalesce(label_facing_alt, label_facing)
  ) %>%
  select(-row_id, -label_y_alt, -label_facing_alt) %>%
  mutate(sector = "plasmid")

circos.clear()

circos.par("gap.degree" = 0, "start.degree" = 1)
circos.initialize(factors = mog$sector, xlim = c(0, 119538))


circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.1, bg.border = "black",
  panel.fun = function(x, y) {
    for (i in 1:nrow(mog)) {
      start <- mog$gene_start[i]
      end <- mog$gene_end[i]
      label <- mog$gene_name[i]
      color <- mog$color[i]
      label_y <- mog$label_y[i]
      label_facing <- mog$label_facing[i]
      
      circos.rect(start, 0, end, 1, col = color, border = "black")
      
      if (!is.na(label) && label != "" && label != "NA") {
        circos.text((start + end) / 2, label_y, labels = label,
                    facing = label_facing, niceFacing = TRUE, cex = 1.2, font = 3)
      }
    }
  }
)

legend("center", legend = names(palette), fill = palette, border = "black", cex = 1.2)

dev.off()

