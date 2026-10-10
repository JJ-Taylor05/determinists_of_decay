# Plot the distribution of half-life values, shading the top and bottom 10%
# ---------------------------------------------------------------------------

library(ggplot2)
library(sysfonts)
library(showtext)

font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)


# ---- 1. Load data ---------------------------------------------------------
# The file has no header: Ensembl ID, gene symbol, half-life value.
# Update the path if the CSV is not in your working directory
# (or use setwd() / Session > Set Working Directory in RStudio).
file_path <- "plot_data/human_halflife_data_sorted.csv"

df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

names(df)[1:3] <- c("ensembl_id", "gene_symbol", "halflife")

df <- df[!is.na(df$halflife), ]

# ---- 2. Work out the 10% / 90% cut-offs -----------------------------------
q_low  <- quantile(df$halflife, 0.10)
q_high <- quantile(df$halflife, 0.90)

df$group <- ifelse(df$halflife <= q_low,  "Bottom 10%",
            ifelse(df$halflife >= q_high, "Top 10%", "Middle 80%"))
df$group <- factor(df$group, levels = c("Bottom 10%", "Middle 80%", "Top 10%"))

# ---- 3. Plot --------------------------------------------------------------
# Histogram with bars coloured by group. Using a fixed set of breaks means no
# bar straddles a cut-off unevenly; the cut-offs are drawn as dashed lines.
p <- ggplot(df, aes(x = halflife, fill = group)) +
  geom_histogram(binwidth = 0.5, boundary = 0, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("Bottom 10%" = "darkorchid4",
                               "Middle 80%" = "grey70",
                               "Top 10%"    = "chocolate3"),
                    name = NULL) +
  labs(title = "Distribution of half-life values",
       subtitle = sprintf("10th percentile = %.2f, 90th percentile = %.2f (n = %d genes)",
                          q_low, q_high, nrow(df)),
       x = "Half-life value",
       y = "Number of genes") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

print(p)

# ---- 4. Save to file -------------------------------------------
ggsave("plots/halflife_distribution.png", p, width = 8, height = 5, dpi = 300)
