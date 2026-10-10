library(tidyverse)
library(Manu)
library(sysfonts)
library(showtext)

# Font & colour setup
font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

kokako <- get_pal("Kokako")

sig_line <- 0.05

# Load in CSVs
summary_df <- read_csv("plot_data/motif_robustness_summary.csv")
null_stable <- read_csv("plot_data/global_null_stable.csv")
null_unstable <- read_csv("plot_data/global_null_unstable.csv")

# Make output directory
output_dir <- "plots/null_control_plots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Plotting functions
# Motif recurrence rate in null, per direction
make_recurrence_plot <- function(df, direction_label, bar_colour) {

  d <- df %>%
    filter(direction == direction_label) %>%
    mutate(
      recurrence_rate = n_null_iterations_with_match / n_iterations_evaluated,
      motif_id = fct_reorder(motif_id, recurrence_rate)
    )

  ggplot(d, aes(x = recurrence_rate, y = motif_id)) +
    geom_col(fill = bar_colour) +
    scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(title = paste0(str_to_title(direction_label), " direction: motif recurrence in null"),
         x = "Raw empirical recurrence rate in null", y = NULL) +
    theme_minimal(base_family = "arimo", base_size = 16) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())
}

stable_recurrence_plot <- make_recurrence_plot(summary_df, "stable", "chocolate3")
stable_recurrence_plot
ggsave(file.path(output_dir, "motif_recurrence_rates_stable.png"), stable_recurrence_plot,
       width = 9, height = 7)

unstable_recurrence_plot <- make_recurrence_plot(summary_df, "unstable", "darkorchid4")
unstable_recurrence_plot
ggsave(file.path(output_dir, "motif_recurrence_rates_unstable.png"), unstable_recurrence_plot,
       width = 9, height = 7)

# Null distribution of motifs passing E-value threshold, per direction
make_null_histogram <- function(df, direction_label) {

  real_n <- summary_df %>% filter(direction == direction_label) %>% nrow()
  pct_extreme <- mean(df$n_motifs <= real_n) * 100

  ggplot(df, aes(x = n_motifs)) +
    geom_histogram(binwidth = 1, fill = kokako[2], colour = "white", boundary = 0) +
    geom_vline(xintercept = real_n, linetype = "dashed", colour = "firebrick", linewidth = 1) +
    annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
             label = paste0("Real run (", real_n, ")\n",
                            sprintf("%.0f%%", pct_extreme), " of null iterations\n",
                            "as/more extreme than real"),
             colour = "firebrick", fill = "white", size = 4) +
    labs(title = paste0(str_to_title(direction_label), "-direction null (n=", nrow(df), ")"),
         x = "Motifs passing E-value threshold", y = "Number of null iterations") +
    theme_minimal(base_family = "arimo", base_size = 16)
}

stable_null_plot <- make_null_histogram(null_stable, "stable")
stable_null_plot
ggsave(file.path(output_dir, "null_distribution_motif_counts_stable.png"), stable_null_plot,
       width = 8, height = 6)

unstable_null_plot <- make_null_histogram(null_unstable, "unstable")
unstable_null_plot
ggsave(file.path(output_dir, "null_distribution_motif_counts_unstable.png"), unstable_null_plot,
       width = 8, height = 6)

