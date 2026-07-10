library(tidyverse)
library(circlize)
library(data.table)
library(ggplotify)

# Figure 2. Plasmid plot
sample_barcode <- "PaPr00000088MP" #change barcode for plasmid generation

mog_results <- fread(file.path("data/plasmid_plot",sample_barcode,"mog.tsv"))
contig_arg <- fread(file.path("data/plasmid_plot",sample_barcode,"contig_arg.tsv"))

# Add degree columns t
plasmid_lengths <- tibble(
  sample_barcode = c("PaPr00000245MP", "PaPr00000001MP", "PaPr00000088MP"),
  plasmid_length = c(90738, 119538, 93918)
)

plasmid_length <- plasmid_lengths %>%
  filter(sample_barcode == !!sample_barcode) %>%
  pull(plasmid_length)



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

plasmid_gg <- as.ggplot(function() {
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
                            #r_outside <- rou2 - 0.06 
                            r_outside <- rou2 - 0.08 

                            for (i in 1:nrow(mog)) {

                                label <- mog$gene_name[i]
                                deg <- mog$deg_end[i] - mog$deg_start[i]

                                if (str_detect(label, paste(c("tra", "rep", "mob"), collapse = "|")) & deg < -4 | label %in% contig_arg$gene_name) {

                                    mid_deg <- (mog$deg_start[i] + mog$deg_end[i]) / 2
                                    mid_rad <- mid_deg * pi / 180

                                    if (label %in% contig_arg$gene_name) {
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
                            y_leg = if (sample_barcode == "PaPr00000088MP") {0.35} else {0.45} #change barcode here
                            legend(x = -0.45, y=y_leg, legend = names(palette), fill = palette, border = "black", cex = 0.8, bg="white", bty = "n")

                            # dev.off()
         })
plasmid_gg


# ggsave(
#   filename = "./plots/PaPr00000088MP_plas.svg",  # or .pdf, .svg, etc.
#   plot = plasmid_gg,
#   width = 13,
#   height = 13,        # adjust based on your needs
#    dpi = 600,          # high resolution for publication
#   units = "cm"
# )
