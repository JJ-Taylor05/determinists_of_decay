library(Manu)
library(readxl)
library(dplyr)
library(ggplot2)
library(sysfonts)
library(showtext)

# Font setup 
font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

# Read and clean data 
corr_data <- read_excel("plot_data/stab_feature_corrs.xlsx") %>%
  filter(species == "human") %>%
  mutate(metric_display = trimws(metric_display))  

# Control display order 
feature_order <- c("Exon junction density", "MFE", "Length", "GC content", "Motifs")
region_order  <- c("5utr", "cds", "3utr", "mrna")
region_labels <- c(`5utr` = "5' UTR", cds = "CDS", `3utr` = "3' UTR", mrna = "mRNA")

corr_data <- corr_data %>%
  mutate(
    metric_display = factor(metric_display, levels = feature_order),
    region = factor(region, levels = region_order, labels = region_labels)
  )

# Alternate panel shading
strip_shade <- data.frame(
  metric_display = factor(feature_order, levels = feature_order),
  shade = rep(c("grey92", "white"), length.out = length(feature_order))
)

# Colour palette 
# First run devtools::install_github("G-Thomson/Manu") if you do not have the native bird colours package
kokako <- get_pal("Kokako")
region_colours <- c(
  "5' UTR" = kokako[5],
  "CDS" = kokako[3],
  "3' UTR" = kokako[2],
  "mRNA" = kokako[1]
)

# Plot
p <- ggplot(corr_data, aes(x = correlation_abs, y = region, colour = region)) +
  geom_rect(
    data = strip_shade,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = shade),
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_errorbarh(aes(xmin = conf.low_abs, xmax = conf.high_abs), height = 0.2) +
  geom_point(size = 2.5) +
  facet_grid(rows = vars(metric_display), scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(values = region_colours) +
  labs(title = "Correlation of mRNA stability features to decay rate", x = "Absolute Spearman correlation", y = NULL, color = "Transcript region") +
  theme_minimal(base_family = "arimo", base_size = 16) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(), 
    plot.title = element_text(hjust = 1.0),
    legend.position = "bottom",
    legend.justification = "left",
    legend.location = "plot"
  )

print(p)

# Save plot
ggsave("stab_feature_corr_plot.png", plot = p, width = 7, height = 6, units = "in", dpi = 300)
