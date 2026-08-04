library(ggplot2)
library(sysfonts)
library(showtext)

font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

stable_path   <- "plot_data/stable_train.fa"
unstable_path <- "plot_data/unstable_train.fa"

read_fasta <- function(path) {
  lines <- readLines(path)
  is_header <- grepl("^>", lines)
  ids <- lines[is_header]
  seq_id <- cumsum(is_header)
  seqs <- vapply(
    split(lines[!is_header], seq_id[!is_header]),
    function(x) paste(x, collapse = ""),
    character(1)
  )
  names(seqs) <- sub("^>", "", ids)
  toupper(seqs)
}

seq_composition <- function(seqs, group_label) {
  bases <- c("A", "C", "G", "T")
  mat <- t(vapply(seqs, function(s) {
    chars <- strsplit(s, "")[[1]]
    total <- length(chars)
    counts <- table(factor(chars, levels = bases))
    as.numeric(counts) / total
  }, numeric(4)))
  colnames(mat) <- bases
  df <- as.data.frame(mat)
  df$seq_id <- names(seqs)
  df$group <- group_label
  df
}

stable_seqs   <- read_fasta(stable_path)
unstable_seqs <- read_fasta(unstable_path)

stable_df   <- seq_composition(stable_seqs,   "Stable")
unstable_df <- seq_composition(unstable_seqs, "Unstable")

comp_df <- rbind(stable_df, unstable_df)
comp_df$group <- factor(comp_df$group, levels = c("Stable", "Unstable"))

bases <- c("A", "C", "G", "T")

test_results <- do.call(rbind, lapply(bases, function(b) {
  x <- comp_df[comp_df$group == "Stable", b]
  y <- comp_df[comp_df$group == "Unstable", b]

  wt <- wilcox.test(x, y)
  tt <- t.test(x, y)

  data.frame(
    nucleotide      = b,
    mean_stable     = mean(x),
    mean_unstable   = mean(y),
    diff            = mean(x) - mean(y),
    wilcox_p        = wt$p.value,
    welch_t_p       = tt$p.value
  )
}))

test_results$wilcox_p_adj  <- p.adjust(test_results$wilcox_p,  method = "BH")
test_results$welch_t_p_adj <- p.adjust(test_results$welch_t_p, method = "BH")

cat("\n---- Per-nucleotide comparison: Stable vs Unstable ----\n")
print(test_results, digits = 4, row.names = FALSE)

sig_stars <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}
test_results$sig <- vapply(test_results$wilcox_p_adj, sig_stars, character(1))

comp_df$GC <- comp_df$C + comp_df$G

gc_test <- wilcox.test(comp_df$GC[comp_df$group == "Stable"],
                        comp_df$GC[comp_df$group == "Unstable"])

gc_sig <- sig_stars(gc_test$p.value)

cat(sprintf("GC%%: Stable mean = %.2f%%, Unstable mean = %.2f%%, Wilcoxon p = %.3g (%s)\n",
            mean(comp_df$GC[comp_df$group == "Stable"]) * 100,
            mean(comp_df$GC[comp_df$group == "Unstable"]) * 100,
            gc_test$p.value, gc_sig))

stack_levels <- c("C", "G", "U", "A") 

plot_df <- do.call(rbind, lapply(bases, function(b) {
  lbl <- if (b == "T") "U" else b
  data.frame(
    group      = c("Stable", "Unstable"),
    nucleotide = lbl,
    pct        = c(mean(comp_df[comp_df$group == "Stable", b]) * 100,
                    mean(comp_df[comp_df$group == "Unstable", b]) * 100)
  )
}))
plot_df$nucleotide <- factor(plot_df$nucleotide, levels = stack_levels)
plot_df$group <- factor(plot_df$group, levels = c("Stable", "Unstable"))

base_colors <- c("A" = "#00811B", "U" = "#CC0000", "G" = "#FF9900", "C" = "#0000CC")

bracket_df <- data.frame(
  y      = 104,
  label  = sprintf("GC %s", gc_sig),
  metric = "GC"
)

p <- ggplot(plot_df, aes(x = group, y = pct, fill = nucleotide)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 2, y = bracket_df$y, yend = bracket_df$y) +
  annotate("segment", x = 1, xend = 1, y = bracket_df$y - 2, yend = bracket_df$y) +
  annotate("segment", x = 2, xend = 2, y = bracket_df$y - 2, yend = bracket_df$y) +
  annotate("text", x = 1.5, y = bracket_df$y + 3, label = bracket_df$label,
           fontface = "bold", size = 5) +
  scale_fill_manual(values = base_colors, breaks = c("A", "U", "G", "C")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 108),
                      breaks = seq(0, 100, 25)) +
  labs(
    title = "Nucleotide composition: Stable vs Unstable transcripts",
    caption = "Wilcoxon rank-sum test, BH-adjusted: *** p<0.001  ** p<0.01  * p<0.05  ns = not significant",
    x = NULL, y = "% of sequence", fill = "Base"
  ) +
  theme_minimal(base_size = 18, base_family = "arimo") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(margin = margin(b = 2)),
    plot.caption  = element_text(size = 12, color = "grey30", hjust = 0)
  )

print(p)

ggsave("nucleotide_freq_comparison.png", p, width = 9, height = 5.5, dpi = 300)

