library(tidyverse)
library(sysfonts)
library(showtext)
library(scales)

font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

region_levels <- c("5' UTR", "CDS", "3' UTR")
region_colours <- c(
  "5' UTR" = "grey50",
  "CDS" = "grey15",
  "3' UTR" = "grey50"
)

stable_path   <- "plot_data/stable_structure_mapping_hits.tsv"
unstable_path <- "plot_data/unstable_structure_mapping_hits.tsv"

prepare_hits <- function(path) {
  read_tsv(path, na = "NA") %>%
    mutate(
      hit_centre = (start + end) / 2,
      norm_position = case_when(
        is.na(orf_start) ~ hit_centre / seq_len,
        hit_centre < orf_start ~ if_else(orf_start > 0, (hit_centre / orf_start) / 3, 0),
        hit_centre < orf_end   ~ 1/3 + if_else((orf_end - orf_start) > 0,
                                               ((hit_centre - orf_start) / (orf_end - orf_start)) / 3, 0),
        TRUE ~ 2/3 + if_else((seq_len - orf_end) > 0,
                             ((hit_centre - orf_end) / (seq_len - orf_end)) / 3, 0)
      )
    )
}

hits_stable   <- prepare_hits(stable_path)
hits_unstable <- prepare_hits(unstable_path)

region_breaks <- c(
  seq(0,   1/3, length.out = 11),
  seq(1/3, 2/3, length.out = 11)[-1],
  seq(2/3, 1,   length.out = 11)[-1]
)
bin_centres <- (head(region_breaks, -1) + tail(region_breaks, -1)) / 2

bin_hits <- function(d, label) {
  bin_index <- cut(d$norm_position, breaks = region_breaks,
                    include.lowest = TRUE, labels = FALSE)
  tibble(bin_index = bin_index) %>%
    filter(!is.na(bin_index)) %>%
    count(bin_index, name = "n_hits") %>%
    complete(bin_index = seq_along(bin_centres), fill = list(n_hits = 0)) %>%
    mutate(
      bin_mid = bin_centres[bin_index],
      proportion = n_hits / sum(n_hits),
      stability = label
    )
}

hits_binned <- bind_rows(
  bin_hits(hits_stable,   "Stabilising motifs"),
  bin_hits(hits_unstable, "Destabilising motifs")
) %>%
  mutate(
    stability = factor(stability, levels = c("Stabilising motifs", "Destabilising motifs")),
    y = if_else(stability == "Destabilising motifs", -proportion, proportion)
  )

stability_colours <- c("Stabilising motifs" = "mediumpurple3", "Destabilising motifs" = "purple4")
y_max <- max(abs(hits_binned$y)) * 1.1

track_gap    <- y_max * 0.15
box_y        <- -(y_max + track_gap * 2.5)
utr_half     <- y_max * 0.06
cds_half     <- y_max * 0.11
box_top      <- box_y + max(utr_half, cds_half)
label_gap    <- track_gap * 0.5
label_y      <- box_top + label_gap
y_lower      <- box_y - cds_half - track_gap

region_shading <- tibble(
  xmin = c(0, 1/3, 2/3), xmax = c(1/3, 2/3, 1),
  region = factor(region_levels, levels = region_levels)
) %>%
  mutate(colour = unname(region_colours[as.character(region)]))

shade_layer <- geom_rect(
  data = region_shading,
  aes(xmin = xmin, xmax = xmax, fill = I(colour)),
  ymin = -y_max, ymax = y_max, alpha = 0.10, inherit.aes = FALSE
)

region_dividers <- tibble(x = c(1/3, 2/3))
divider_layer <- geom_segment(
  data = region_dividers,
  aes(x = x, xend = x, y = y_max, yend = box_top),
  colour = "grey70", linetype = "dashed", inherit.aes = FALSE
)

transcript_cartoon <- tibble(
  region = factor(region_levels, levels = region_levels),
  xmin = c(0, 1/3, 2/3),
  xmax = c(1/3, 2/3, 1),
  half_height = c(utr_half, cds_half, utr_half)
) %>%
  mutate(colour = unname(region_colours[as.character(region)]))

transcript_layer <- geom_rect(
  data = transcript_cartoon,
  aes(xmin = xmin, xmax = xmax, ymin = box_y - half_height, ymax = box_y + half_height, fill = I(colour)),
  colour = "grey10", linewidth = 0.3,
  inherit.aes = FALSE
)

region_label_layer <- annotate(
  "text", x = c(1/6, 0.5, 5/6), y = rep(label_y, 3),
  label = region_levels, family = "arimo", size = 5, colour = "grey30"
)

p_metatranscript <- ggplot(hits_binned, aes(x = bin_mid, y = y, colour = stability, fill = stability)) +
  shade_layer + divider_layer +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_area(alpha = 0.35, position = "identity") +
  geom_line(linewidth = 0.9) +
  region_label_layer +
  transcript_layer +
  annotate("text", x = 1/6, y = y_max * 0.55, label = "Stable-associated structure",
           colour = stability_colours[["Stabilising motifs"]],
           family = "arimo", size = 4) +
  annotate("text", x = 1/6, y = -y_max * 0.55, label = "Unstable-associated structure",
           colour = stability_colours[["Destabilising motifs"]],
           family = "arimo", size = 4) +
  scale_colour_manual(values = stability_colours, name = NULL) +
  scale_fill_manual(values = stability_colours, name = NULL) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(y_lower, y_max), expand = expansion(mult = c(0, 0.02)),
    breaks = function(x) { b <- scales::extended_breaks()(c(-y_max, y_max)); b[abs(b) <= y_max] },
    labels = function(v) number(abs(v), accuracy = 0.01)
  ) +
  labs(
    x = "Meta-transcript position",
    y = "Proportion of transcripts with structure",
    title = "Meta-transcript distribution of structured regions"
  ) +
  theme_minimal(base_family = "arimo", base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(hjust = (0 - y_lower) / (y_max - y_lower))
  )

p_metatranscript

ggsave(("metatranscript_stable_vs_unstable_structure.png"), p_metatranscript,
       width = 10, height = 6, dpi = 300, bg = "white")

