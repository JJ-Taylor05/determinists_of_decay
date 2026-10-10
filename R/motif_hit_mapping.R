# R version of the python script that maps motif hits
# back to regions of the transcript

library(tidyverse)
library(Manu)
library(sysfonts)
library(showtext)
library(patchwork)
library(scales)

font_add_google("Arimo", "arial")
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


set.seed(1)      # makes the random control reproducible
n_sim <- 5000    # number of random "worlds" to simulate

# Helper functions
hits <- hits %>%
  mutate(
    motif_num = as.integer(factor(motif_label)),               # motif -> 1..n
    group_num = as.integer(factor(paste(seq_id, motif_id)))    # one id per (transcript, motif) pair
  )
motif_names <- levels(factor(hits$motif_label))
n_motifs    <- length(motif_names)
n_groups    <- max(hits$group_num)
n_slots     <- n_motifs * 3     # one slot per motif x region

# Random simulation
null_counts <- matrix(0L, nrow = n_sim, ncol = n_slots)

for (i in seq_len(n_sim)) {
  shift      <- runif(n_groups)[hits$group_num] * hits$seq_len        # one random shift per (transcript, motif)
  new_centre <- (hits$hit_centre + shift) %% hits$seq_len             # slide hits round the transcript
  new_region <- 1L + (new_centre >= hits$orf_start) + (new_centre >= hits$orf_end)  # 1 = 5' UTR, 2 = CDS, 3 = 3' UTR
  null_counts[i, ] <- tabulate((hits$motif_num - 1L) * 3L + new_region, nbins = n_slots)
}

# Random summarize
null_summary <- tibble(
  motif_label = rep(motif_names, each = 3),
  region      = factor(rep(region_levels, times = n_motifs), levels = region_levels),
  exp_mean    = colMeans(null_counts),
  exp_lo      = apply(null_counts, 2, quantile, probs = 0.025),
  exp_hi      = apply(null_counts, 2, quantile, probs = 0.975)
) %>%
  left_join(region_counts, by = c("motif_label", "region")) %>%
  mutate(
    p_high  = (colSums(sweep(null_counts, 2, n_hits, ">=")) + 1) / (n_sim + 1),
    p_low   = (colSums(sweep(null_counts, 2, n_hits, "<=")) + 1) / (n_sim + 1),
    p_value = pmin(1, 2 * pmin(p_high, p_low)),
    p_adj   = p.adjust(p_value, method = "BH"),
    log2_obs_vs_exp = log2((n_hits + 0.5) / (exp_mean + 0.5))
  )

# Plot
p_counts <- ggplot(region_counts, aes(x = motif_label, y = n_hits, fill = region)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_linerange(data = null_summary,
                 aes(y = exp_mean, ymin = pmax(exp_lo, 0.5), ymax = exp_hi),
                 position = position_dodge(width = 0.7), colour = "black", linewidth = 0.4) +
  geom_errorbar(data = null_summary,
                aes(y = exp_mean, ymin = exp_mean, ymax = exp_mean),
                position = position_dodge(width = 0.7), width = 0.65,
                colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = region_colours, name = "Region") +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Motifs enriched in unstable mRNA",
       y = expression("Number of hits (log"[10]*")"),
       title = "Hit counts per region for motifs enriched in unstable mRNA",
       caption = "Black tick = mean expected hits under random placement; black line = 95% interval (5,000 simulations)") +
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
  labs(x = NULL, y = "Hits (%)", title = "Regional distribution of unstable motif hits") +
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
    theme_minimal(base_family = "arimo", base_size = 30)
  
  n_label <- ggplot() +
    theme_void() +
    annotate(
      "text", x = 0.9, y = 0.9,
      label = str_wrap(paste0("Total hits = ", nrow(d)), width = 12),
      family = "arimo", size = 10, lineheight = 0.9,
      hjust = 1, vjust = 1,
      colour = "grey20"
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1))
  
  p_hist + inset_element(
    n_label,
    left = 0.78, bottom = 0.7, right = 1.0, top = 1.0,
    align_to = "full"
  )
}

motif_ids <- sort(unique(hits$motif_id))
motif_pages <- map(motif_ids, make_motif_page)

motif_pages[1]

# Create output directory
out_dir <- "plots/unstable_motifmapping_plots"
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
