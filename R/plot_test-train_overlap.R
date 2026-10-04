## Reads: test_vs_train.tsv  - from mmseqs easy-search (test as query, train as target)
## Set working directory to test_train_similarity_check/

library(ggplot2)

# Change depending on stable/unstable set
search_tsv <- "test_vs_train_unstable.tsv"

# ---- motif-discovery-aware leakage thresholds ----
# Tuned for motifs in the 10-20nt range. A real motif can legitimately recur
# in both train and test at high identity over ~10-20bp - that's the whole
# point of a generalizable motif, not leakage. What IS leakage is a much
# longer near-identical stretch, since that means the *same sequence
# instance* (not just the same motif) is sitting in both sets.
longest_motif_nt   <- 20
min_alnlen_partial <- max(50, 5 * longest_motif_nt)   # 100bp here
min_pident_partial <- 90   # identity required for the length-based check
min_pident_whole    <- 80   # identity required for the whole-duplicate check
min_cov_whole       <- 0.8  # coverage required on AT LEAST ONE side (not both -
                             # requiring both misses a shorter isoform fully
                             # contained in a longer transcript)
# ---------------------------------------------------

hit_cols <- c("query", "target", "pident", "alnlen", "mismatch", "gapopen",
              "qstart", "qend", "tstart", "tend", "evalue", "bits",
              "qlen", "tlen", "qcov", "tcov")

hits <- read.delim(search_tsv, header = FALSE, col.names = hit_cols,
                    stringsAsFactors = FALSE)

if (nrow(hits) == 0) {
  stop("No hits found in ", search_tsv, " - nothing to plot.")
}

hits$is_whole_dup <- hits$pident >= min_pident_whole &
  (hits$qcov >= min_cov_whole | hits$tcov >= min_cov_whole)
hits$is_partial_dup <- hits$pident >= min_pident_partial &
  hits$alnlen >= min_alnlen_partial

hits$flag <- "none"
hits$flag[hits$is_partial_dup] <- "partial_duplicate"
hits$flag[hits$is_whole_dup]   <- "whole_duplicate"  # whole takes priority if both trip
hits$flag <- factor(hits$flag, levels = c("none", "partial_duplicate", "whole_duplicate"))

flagged <- hits[hits$flag != "none", ]
flagged <- flagged[order(-flagged$alnlen), ]
if (nrow(flagged) > 0) {
  write.csv(flagged, "leakage_candidates.csv", row.names = FALSE)
  cat(sprintf("%d hit(s) flagged as likely leakage (see leakage_candidates.csv)\n", nrow(flagged)))
  cat(sprintf("  whole_duplicate   (id>=%d%%, cov>=%d%% on either side): %d\n",
              min_pident_whole, min_cov_whole * 100, sum(flagged$flag == "whole_duplicate")))
  cat(sprintf("  partial_duplicate (id>=%d%%, alnlen>=%dbp, any coverage): %d\n",
              min_pident_partial, min_alnlen_partial, sum(flagged$flag == "partial_duplicate")))
} else {
  cat("No hits flagged as likely leakage at the current thresholds.\n")
}

hits$cov_either <- pmax(hits$qcov, hits$tcov) * 100  # matches the whole_duplicate rule exactly

pB <- ggplot(hits, aes(x = cov_either, y = pident, color = flag)) +
  geom_jitter(width = 0.5, height = 0.5, alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c(none = "#8C9196", partial_duplicate = "#DD8452", whole_duplicate = "#C44E52"),
    labels = c(none = "below threshold", partial_duplicate = "partial duplicate",
               whole_duplicate = "whole duplicate"),
    name = NULL, drop = FALSE) +
  coord_cartesian(xlim = c(0, 100), ylim = c(min(60, min(hits$pident) - 2), 102)) +
  labs(title = "Identity vs. coverage of test\u2192train hits",
       x = "Coverage (%)", y = "Percent identity (%)") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20),
        legend.position = "bottom", legend.text = element_text(size = 12))

pB

ggsave("identity_vs_coverage.png", pB, width = 7, height = 6.5, dpi = 150)
