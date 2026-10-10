# ── Packages ──────────────────────────────────────────────────────────────────
library(ggplot2)
library(ggrepel)

# ── Data ──────────────────────────────────────────────────────────────────────
df <- data.frame(
  tool      = c("MEME", "HOMER", "ChIPMunk", "Autoseed", "STREME", "Dimont", "ExplaiNN", "RCade", "gkmSVM", "ProBound"),   # Replace with your tool names
  rank      = c(4, 6, 2, 5, 7, 1, 8, 9, 10, 3),                         # Replace with your rank scores (1–10)
  citations = c(327, 843, 9.9, 36.4, 55, 3.5, 33.3, 7.3, 24.8, 27.8)                     # Replace with your citation counts
)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(df, aes(x = rank, y = citations)) +
  
  # Scatter points coloured by citation count
  geom_point(aes(colour = rank), size = 4) +
  
  # Tool name labels that repel each other to avoid overlap
  geom_text_repel(
    aes(label = tool),
    colour         = "black",
    size           = 3.5,
    fontface       = "italic",
    box.padding    = 0.8,
    point.padding  = 0.3,
    max.overlaps   = Inf,
    force          = 2,
    segment.colour = "transparent"
  ) +
  
  # Colour gradient for points and legend
  scale_colour_gradient(
    low  = "#06948E",
    high = "#C84D4C",
    name = "Citations"
  ) +
  
  # Ensure all rank values 1–10 appear on the x axis
  scale_y_log10(
    breaks = 10^(0:3),
    labels = function(x) parse(text = paste0("10^", round(log10(x))))
  ) +
  scale_x_continuous(breaks = 1:10) +
  
  # Axis and title labels
  labs(
    title = "Motif Finding Tools: Rank vs Popularity",
    x     = "Relative benchmark score",
    y     = "Citations per year"
  ) +
  
  # Clean theme with font styling
  theme_classic() +
  theme(
    plot.title   = element_text(size = 16, face = "bold",   hjust = 0.5),
    axis.title   = element_text(size = 13, face = "italic"),
    axis.text    = element_text(size = 11),
    legend.title = element_text(face = "bold")
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("plots/motif_tools_plot.png", width = 8, height = 4, dpi = 300)
