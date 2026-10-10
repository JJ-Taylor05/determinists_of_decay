library(ggplot2)
library(showtext)

font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

base_family <- "arimo"

# Read data
stable_path   <- "plot_data/stable_train.fa"
unstable_path <- "plot_data/unstable_train.fa"
merged_path   <- "plot_data/merged_scores.tsv"  

group_colors <- c("Stable" = "chocolate3", "Unstable" = "darkorchid4")

# Shared legend key: a plain coloured square with a black outline
draw_key_outlined_square <- function(data, params, size) {
  grid::rectGrob(
    width  = grid::unit(0.8, "npc"),
    height = grid::unit(0.8, "npc"),
    gp = grid::gpar(col = "black",
                    fill = scales::alpha(data$fill, 0.85),
                    lwd = 1, linejoin = "mitre")
  )
}

read_fasta <- function(path) {
  lines     <- readLines(path)
  is_header <- grepl("^>", lines)
  ids       <- lines[is_header]
  seq_id    <- cumsum(is_header)
  seqs <- vapply(
    split(lines[!is_header], seq_id[!is_header]),
    function(x) paste(x, collapse = ""),
    character(1)
  )
  names(seqs) <- trimws(sub("^>", "", ids))
  toupper(seqs)
}

## GC, bitscore and half-life comparison
# Seq composition function
seq_composition <- function(seqs, group_label) {
  bases <- c("A", "C", "G", "T")
  mat <- t(vapply(seqs, function(s) {
    chars  <- strsplit(s, "")[[1]]
    counts <- table(factor(chars, levels = bases))
    as.numeric(counts) / length(chars)
  }, numeric(4)))
  colnames(mat) <- bases
  df <- as.data.frame(mat)
  df$transcript_id <- names(seqs)
  df$group <- group_label
  df
}

sig_stars <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

## Nucleotide composition from the FASTA files
comp_df <- rbind(
  seq_composition(read_fasta(stable_path),   "Stable"),
  seq_composition(read_fasta(unstable_path), "Unstable")
)
comp_df$group <- factor(comp_df$group, levels = c("Stable", "Unstable"))
comp_df$GC    <- comp_df$C + comp_df$G

## Half-life and bitscore sums, joined by transcript ID
scores <- read.delim(merged_path, stringsAsFactors = FALSE, strip.white = TRUE)
scores <- scores[, c("transcript_id", "degradation_value", "sum_of_bits_combined_score")]

all_df <- merge(comp_df, scores, by = "transcript_id", all.x = TRUE)

## Sanity check: how many sequences have score data?
message(sprintf("Sequences: %d stable, %d unstable",
                sum(all_df$group == "Stable"), sum(all_df$group == "Unstable")))
message(sprintf("With half-life/bitscore data: %d stable, %d unstable",
                sum(all_df$group == "Stable"   & !is.na(all_df$degradation_value)),
                sum(all_df$group == "Unstable" & !is.na(all_df$degradation_value))))


# Reformat
metrics <- list(
  GC       = list(col = "GC",                         scale = 100, label = "GC (%)"),
  HalfLife = list(col = "degradation_value",          scale = 1,   label = "Half-life value"),
  Bitscore = list(col = "sum_of_bits_combined_score", scale = 1,   label = "Motif bitscore sum")
)

long_df <- do.call(rbind, lapply(names(metrics), function(m) {
  d <- data.frame(
    group  = all_df$group,
    metric = m,
    value  = all_df[[metrics[[m]]$col]] * metrics[[m]]$scale
  )
  d[!is.na(d$value), ]
}))
long_df$metric <- factor(
  long_df$metric,
  levels = names(metrics),
  labels = vapply(metrics, `[[`, character(1), "label")
)


# Statistics
## Wilcoxon rank-sum test per metric, BH-adjusted across all panels.
## rank_biserial: effect size from -1 to 1 (0 = no separation, +/-1 = complete).
stats_df <- do.call(rbind, lapply(levels(long_df$metric), function(m) {
  d  <- long_df[long_df$metric == m, ]
  x  <- d$value[d$group == "Stable"]
  y  <- d$value[d$group == "Unstable"]
  wt <- wilcox.test(x, y)
  data.frame(
    metric          = m,
    n_stable        = length(x),
    n_unstable      = length(y),
    median_stable   = median(x),
    median_unstable = median(y),
    rank_biserial   = 2 * unname(wt$statistic) / (length(x) * length(y)) - 1,
    wilcox_p        = wt$p.value,
    y_max           = max(d$value),
    y_min           = min(d$value)
  )
}))
stats_df$wilcox_p_adj <- p.adjust(stats_df$wilcox_p, method = "BH")
stats_df$sig          <- vapply(stats_df$wilcox_p_adj, sig_stars, character(1))
stats_df$metric       <- factor(stats_df$metric, levels = levels(long_df$metric))

## Positions for the significance brackets (one per panel)
stats_df$rng    <- stats_df$y_max - stats_df$y_min
stats_df$y_br   <- stats_df$y_max + 0.06 * stats_df$rng
stats_df$y_tick <- stats_df$y_br  - 0.02 * stats_df$rng
stats_df$y_txt  <- stats_df$y_br  + 0.01 * stats_df$rng

## Results table 
results_table <- stats_df[, c("metric", "n_stable", "n_unstable", "median_stable",
                              "median_unstable", "rank_biserial",
                              "wilcox_p", "wilcox_p_adj", "sig")]
View(results_table)


# Plot 
motifclass_comparison <- ggplot(long_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5, outlier.alpha = 0.35,
               linewidth = 0.5, alpha = 0.85, key_glyph = draw_key_outlined_square) +
  geom_segment(data = stats_df, inherit.aes = FALSE,
               aes(x = 1, xend = 2, y = y_br,   yend = y_br)) +
  geom_segment(data = stats_df, inherit.aes = FALSE,
               aes(x = 1, xend = 1, y = y_tick, yend = y_br)) +
  geom_segment(data = stats_df, inherit.aes = FALSE,
               aes(x = 2, xend = 2, y = y_tick, yend = y_br)) +
  geom_text(data = stats_df, inherit.aes = FALSE,
            aes(x = 1.5, y = y_txt, label = sig),
            vjust = 0, fontface = "bold", size = 5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  labs(
    title   = "Stable vs Unstable transcripts: GC content, half-life and motif scores",
    caption = paste0("Wilcoxon rank-sum test, BH-adjusted across all panels: ",
                     "*** p<0.001  ** p<0.01  * p<0.05  ns = not significant\n",
                     "Boxes: median and IQR; whiskers: 1.5 x IQR; points: outliers"),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 16, base_family = base_family) +
  theme(
    legend.position    = "right",
    legend.text        = element_text(size = 14),
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text    = element_text(face = "bold", size = 15),
    plot.title    = element_text(margin = margin(b = 6)),
    plot.caption  = element_text(size = 11, color = "grey30", hjust = 0),
    panel.spacing = unit(1.1, "lines")
  )

motifclass_comparison  


# Save 
ggsave("plots/motifclass_comparison.png", motifclass_comparison, width = 11, height = 4.5, dpi = 300, bg = "white")


## Length comparison
boundaries_path <- "plot_data/transcript_boundaries.csv"

# Strip the Ensembl version suffix 
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

# Sequence lengths straight from the FASTA files
stable_seqs   <- read_fasta(stable_path)
unstable_seqs <- read_fasta(unstable_path)

fasta_len <- rbind(
  data.frame(transcript_id = names(stable_seqs),
             fasta_len = nchar(stable_seqs),   group = "Stable"),
  data.frame(transcript_id = names(unstable_seqs),
             fasta_len = nchar(unstable_seqs), group = "Unstable")
)
fasta_len$sequence_id <- strip_version(fasta_len$transcript_id)
fasta_len$group <- factor(fasta_len$group, levels = c("Stable", "Unstable"))

# Region boundaries, joined by (unversioned) transcript ID
bounds <- read.csv(boundaries_path, stringsAsFactors = FALSE)
bounds$has_orf <- as.logical(bounds$has_orf)

len_df <- merge(fasta_len, bounds, by = "sequence_id")
n_matched <- table(len_df$group)

# Keep only transcripts whose boundary-file length equals the FASTA length;
# a mismatch suggests the boundaries came from a different transcript version
len_df <- len_df[len_df$fasta_len == len_df$seq_len, ]

# Transcripts with no predicted ORF have no real UTR/CDS split, so blank
# those three regions (their total length is still used)
no_orf <- !len_df$has_orf
len_df$utr5_len[no_orf] <- NA
len_df$cds_len[no_orf]  <- NA
len_df$utr3_len[no_orf] <- NA

message(sprintf("Matched to boundaries file: %d stable, %d unstable (of %d, %d)",
                n_matched["Stable"], n_matched["Unstable"],
                sum(fasta_len$group == "Stable"), sum(fasta_len$group == "Unstable")))
message(sprintf("After length check: %d stable, %d unstable",
                sum(len_df$group == "Stable"), sum(len_df$group == "Unstable")))
message(sprintf("Without a predicted ORF (excluded from UTR/CDS panels): %d stable, %d unstable",
                sum(no_orf & len_df$group == "Stable"),
                sum(no_orf & len_df$group == "Unstable")))

# Long format: one row per transcript per length measure
len_metrics <- list(
  Total = list(col = "seq_len",  label = "mRNA"),
  UTR5  = list(col = "utr5_len", label = "5' UTR"),
  CDS   = list(col = "cds_len",  label = "CDS"),
  UTR3  = list(col = "utr3_len", label = "3' UTR")
)

long_len <- do.call(rbind, lapply(names(len_metrics), function(m) {
  d <- data.frame(group  = len_df$group,
                  metric = m,
                  value  = len_df[[len_metrics[[m]]$col]])
  d[!is.na(d$value), ]
}))
long_len$metric <- factor(long_len$metric,
                          levels = names(len_metrics),
                          labels = vapply(len_metrics, `[[`, character(1), "label"))

# Statistics: Wilcoxon rank-sum per measure, BH-adjusted across these 4 panels
len_stats <- do.call(rbind, lapply(levels(long_len$metric), function(m) {
  d  <- long_len[long_len$metric == m, ]
  x  <- d$value[d$group == "Stable"]
  y  <- d$value[d$group == "Unstable"]
  wt <- wilcox.test(x, y)
  data.frame(
    metric          = m,
    n_stable        = length(x),
    n_unstable      = length(y),
    median_stable   = median(x),
    median_unstable = median(y),
    rank_biserial   = 2 * unname(wt$statistic) / (length(x) * length(y)) - 1,
    wilcox_p        = wt$p.value
  )
}))
len_stats$wilcox_p_adj <- p.adjust(len_stats$wilcox_p, method = "BH")
len_stats$sig          <- vapply(len_stats$wilcox_p_adj, sig_stars, character(1))
len_stats$metric       <- factor(len_stats$metric, levels = levels(long_len$metric))

# The plot uses a log10 axis (lengths are strongly right-skewed), so
# zero-length regions cannot be drawn; they are still counted in the tests
n_zero <- sum(long_len$value <= 0)
if (n_zero > 0) message(n_zero, " zero-length values not shown on the log axis")
plot_len <- long_len[long_len$value > 0, ]

# Bracket positions, computed on the log10 scale
pos <- do.call(rbind, lapply(levels(plot_len$metric), function(m) {
  v <- log10(plot_len$value[plot_len$metric == m])
  data.frame(metric = m, lmax = max(v), lrng = max(v) - min(v))
}))
len_stats <- merge(len_stats, pos, by = "metric")
len_stats$metric <- factor(len_stats$metric, levels = levels(long_len$metric))
len_stats$y_br   <- 10^(len_stats$lmax + 0.06 * len_stats$lrng)
len_stats$y_tick <- 10^(len_stats$lmax + 0.04 * len_stats$lrng)
len_stats$y_txt  <- 10^(len_stats$lmax + 0.07 * len_stats$lrng)

length_table <- len_stats[order(len_stats$metric),
                          c("metric", "n_stable", "n_unstable", "median_stable",
                            "median_unstable", "rank_biserial",
                            "wilcox_p", "wilcox_p_adj", "sig")]
View(length_table)

# Plot
length_comparison <- ggplot(plot_len, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5, outlier.alpha = 0.35,
               linewidth = 0.5, alpha = 0.85, key_glyph = draw_key_outlined_square) +
  geom_segment(data = len_stats, inherit.aes = FALSE,
               aes(x = 1, xend = 2, y = y_br,   yend = y_br)) +
  geom_segment(data = len_stats, inherit.aes = FALSE,
               aes(x = 1, xend = 1, y = y_tick, yend = y_br)) +
  geom_segment(data = len_stats, inherit.aes = FALSE,
               aes(x = 2, xend = 2, y = y_tick, yend = y_br)) +
  geom_text(data = len_stats, inherit.aes = FALSE,
            aes(x = 1.5, y = y_txt, label = sig),
            vjust = 0, fontface = "bold", size = 5) +
  facet_wrap(~ metric, ncol = 4) +   # shared y-axis so regions are directly comparable
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_log10(expand = expansion(mult = c(0.04, 0.14)),
                breaks = 10^(0:5),
                labels = scales::label_math(10^.x, format = log10)) +
  labs(
    title   = "Stable vs Unstable transcripts: length comparisons",
    caption = paste0("Wilcoxon rank-sum test, BH-adjusted across panels: ",
                     "*** p<0.001  ** p<0.01  * p<0.05  ns = not significant\n",
                     "Shared y-axis on a log scale. Boxes: median and IQR; whiskers: 1.5 x IQR; points: outliers"),
    x = NULL, y = expression("Length (nt)")
  ) +
  theme_minimal(base_size = 16, base_family = base_family) +
  theme(
    legend.position    = "right",
    legend.text        = element_text(size = 14),
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text    = element_text(face = "bold", size = 15),
    plot.title    = element_text(margin = margin(b = 6)),
    plot.caption  = element_text(size = 11, color = "grey30", hjust = 0),
    panel.spacing = unit(1.1, "lines")
  )

length_comparison

ggsave("plots/motifclass_length.png", length_comparison, width = 13, height = 4.8, dpi = 300, bg = "white")
