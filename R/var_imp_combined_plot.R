library(tidyverse)
library(Manu)
library(sysfonts)
library(showtext)

# Font & colour setup (matches randomforest.R)
font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

kokako <- get_pal("Kokako")

# --- Config: each file's full path and its region label ---------------------
# Update these paths to wherever each CSV actually lives.
paths <- c(
  "rf_results/all/all_motif_weights.csv",
  "rf_results/cds/cds_motif_weights.csv",
  "rf_results/utr5/utr5_motif_weights.csv",
  "rf_results/utr3/utr3_motif_weights.csv"
)

regions <- c(
  "mRNA",
  "CDS",
  "5'UTR",
  "3'UTR"
)

# --- Load and combine all 4 files -------------------------------------------
combined <- map2_dfr(paths, regions, function(path, region) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(region = region)
})

# --- Build motif-region label ------------------------------------------------
# motif_label in the CSVs stores the sign as a suffix (e.g. "YCAAUAAA+").
# Move the sign to the front and append the region, e.g. "+YCAAUAAA-3'UTR"
combined <- combined %>%
  mutate(
    sign = if_else(str_sub(motif_label, -1) %in% c("+", "-"),
                   str_sub(motif_label, -1), ""),
    motif_stripped = if_else(sign != "",
                             str_sub(motif_label, 1, -2),
                             motif_label),
    plot_label = paste0(sign, motif_stripped, "-", region)
  ) %>%
  dplyr::select(-sign, -motif_stripped)

# --- Select top 20 across all four files -------------------------------------
top20 <- combined %>%
  arrange(desc(importance)) %>%
  slice_head(n = 20)

# --- Plot (styled to match var_imp_plot in randomforest.R) -------------------
top20_plot <- ggplot(top20, aes(x = reorder(plot_label, importance), y = importance)) +
  geom_col(fill = kokako[3]) +
  coord_flip() +
  labs(title = "Top 20 Motif Importances Across All Regions",
       x = "Motifs", y = "Importance") +
  theme_minimal(base_family = "arimo", base_size = 16) +
  theme(plot.title = element_text(hjust = 2.0))

top20_plot
ggsave("top20_motif_importance_plot.png", top20_plot, width = 7, height = 6)
