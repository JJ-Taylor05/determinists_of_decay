# R version of the python script that maps motif hits
# back to regions of the transcript

library(tidyverse)
library(Manu)
library(sysfonts)
library(showtext)
library(patchwork)
library(scales)

font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

kokako <- get_pal("Kokako")

region_levels <- c("5' UTR", "CDS", "3' UTR")
region_colours <- c(
  "5' UTR" = kokako[5],
  "CDS" = kokako[3],
  "3' UTR" = kokako[2]
)

hits_path <- "plot_data/unstable_motif_mapping_hits.tsv"

hits <- read_tsv(hits_path, na = "NA") %>%
  mutate(
    motif_label = str_remove(motif_id, "^\\d+-"),
    hit_centre = (start + end) / 2,
    region = case_when(
      is.na(orf_start)       ~ NA_character_,
      hit_centre < orf_start ~ "5' UTR",
      hit_centre < orf_end   ~ "CDS",
      TRUE                   ~ "3' UTR"
    ),
    region = factor(region, levels = region_levels),
    norm_position = case_when(
      is.na(orf_start) ~ hit_centre / seq_len,
      hit_centre < orf_start ~ if_else(orf_start > 0, (hit_centre / orf_start) / 3, 0),
      hit_centre < orf_end   ~ 1/3 + if_else((orf_end - orf_start) > 0,
                                             ((hit_centre - orf_start) / (orf_end - orf_start)) / 3, 0),
      TRUE ~ 2/3 + if_else((seq_len - orf_end) > 0,
                           ((hit_centre - orf_end) / (seq_len - orf_end)) / 3, 0)
    )
  )

# Grouped bar of raw hit counts per motif per region
region_counts <- hits %>%
  count(motif_label, region, name = "n_hits") %>%
  complete(motif_label, region, fill = list(n_hits = 0))

p_counts <- ggplot(region_counts, aes(x = motif_label, y = n_hits, fill = region)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = region_colours, name = "Region") +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Motifs enriched in unstable mRNA", y = "Number of hits (log10)", title = "Hit counts per region for motifs enriched in unstable mRNA") +
  theme_minimal(base_family = "arimo", base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

p_counts

# Stacked bar chart
region_pct <- region_counts %>%
  group_by(motif_label) %>%
  mutate(pct = n_hits / sum(n_hits) * 100) %>%
  ungroup()

p_pct <- ggplot(region_pct, aes(x = motif_label, y = pct, fill = region)) +
  geom_col(width = 0.55) +
  scale_fill_manual(values = region_colours, name = "Region") +
  labs(x = NULL, y = "Hits (%)", title = "Regional distribution of stable motif hits") +
  theme_minimal(base_family = "arimo", base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 1.0)
  )

p_pct

# Per motif position density
region_shading <- tibble(
  xmin = c(0, 1/3, 2/3), xmax = c(1/3, 2/3, 1),
  region = factor(c("5' UTR", "CDS", "3' UTR"), levels = region_levels)
)
shade_layer <- geom_rect(data = region_shading, aes(xmin = xmin, xmax = xmax, fill = region),
                         ymin = -Inf, ymax = Inf, alpha = 0.10, inherit.aes = FALSE)
divider_layer <- geom_vline(xintercept = c(1/3, 2/3), colour = "grey70", linetype = "dashed")

region_breaks <- c(
  seq(0,   1/3, length.out = 11),
  seq(1/3, 2/3, length.out = 11)[-1],
  seq(2/3, 1,   length.out = 11)[-1]
)

make_motif_page <- function(mid) {
  d <- filter(hits, motif_id == mid)
  lbl <- unique(d$motif_label)
  
  p_hist <- ggplot(d, aes(x = norm_position, fill = region)) +
    shade_layer + divider_layer +
    geom_histogram(aes(y = after_stat(count / sum(count))),
                   breaks = region_breaks, colour = "white", linewidth = 0.2, alpha = 0.9) +
    scale_fill_manual(values = region_colours, name = "Region") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Meta-transcript position", y = "Proportion of hits", title = lbl) +
    theme_minimal(base_family = "arimo", base_size = 25)
  
  n_label <- ggplot() +
    theme_void() +
    annotate(
      "text", x = 0, y = 0,
      label = paste0("Total hits = ", nrow(d), " hits"),
      family = "arimo", size = 5.5, lineheight = 0.9,
      colour = "grey20"
    )
  
  p_hist + inset_element(
    n_label,
    left = 0.82, bottom = 0.90, right = 1.0, top = 1.0,
    align_to = "full"
  )
}

motif_ids <- sort(unique(hits$motif_id))
motif_pages <- map(motif_ids, make_motif_page)

motif_pages[1]

# Create output directory
out_dir <- "unstable_motifmapping_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save the two summary plots
ggsave(file.path(out_dir, "hit_counts_per_region.png"), p_counts,
       width = 14, height = 7, dpi = 300, bg = "white")

ggsave(file.path(out_dir, "hit_pct_per_region.png"), p_pct,
       width = 8, height = 5, dpi = 300, bg = "white")

# Save each per-motif page, named by motif_id
walk2(motif_pages, motif_ids, function(p, mid) {
  ggsave(
    filename = file.path(out_dir, paste0("motif_", mid, "_position.png")),
    plot = p,
    width = 12, height = 6, dpi = 300, bg = "white"
  )
})
