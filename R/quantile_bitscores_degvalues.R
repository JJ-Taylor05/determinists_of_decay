library(tidyverse)
library(boot)
library(Manu)
library(sysfonts)
library(showtext)

# Font & colour setup 
font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

kokako <- get_pal("Kokako")  

## To plot the correlation between deg values and motifs when different filters are applied, simply set your 
## working directory to the folder containing the appropriate 5 quantile files to load in, then adjust the
## labels on the plots.

# Read in quantile files - change these to match filters used
quantile_files = c("q1" = "q1_utr5_merged.tsv",
                   "q2" = "q2_utr5_merged.tsv",
                   "q3" = "q3_utr5_merged.tsv",
                   "q4" = "q4_utr5_merged.tsv",
                   "q5" = "q5_utr5_merged.tsv")

# Putting all data into one dataframe with a column for quantile name
all_data = imap_dfr(quantile_files, function(path, qname) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(quantile = qname)
})

# Calc Spearman correlation for quantile
spearman_stat = function(data, indices) {
  d <- data[indices, ]
  cor(d$degradation_value, d$sum_of_bits_combined_score, method = "spearman")
}

# Make random resampling reproducible
set.seed(123)

# Bootstrapping & calc spearman + CI for all data
boot_all = boot(data = all_data, statistic = spearman_stat, R = 2000)
ci_overall = boot.ci(boot_all, type = "perc")
    
overall_result = tibble(
  rho = boot_all$t0,
  ci_lower = ci_overall$percent[4],
  ci_upper = ci_overall$percent[5],
  n = nrow(all_data)
)

# Scatter plot
scatter_plot = ggplot(all_data, aes(x = degradation_value, y = sum_of_bits_combined_score)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, colour = kokako[2]) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = sprintf("rho = %.2f\n[%.2f, %.2f]\nn = %d", 
                    overall_result$rho, 
                    overall_result$ci_lower, 
                    overall_result$ci_upper,
                    overall_result$n),
    hjust = -0.1, vjust = 1.2, size = 3.0
  ) +
  # Change labels to match filters used
  labs(
    x = "Degradation value",
    y = "Sum of bits",
    colour = "Quantile",
    title = "Degradation values vs sum of bits score in 5' UTR"
  ) +
  theme_bw(base_family = "arimo")

# Quantile level boxplot
all_data = all_data %>%
  mutate(quantile = factor(quantile, levels = c("q1", "q2", "q3", "q4", "q5")))

# Statistical test for trend across ordered groups (Kruskal-Wallis test)
kw_test = kruskal.test(sum_of_bits_combined_score ~ quantile, data = all_data)
print(kw_test)

kw_label = sprintf(
  "Kruskal-Wallis chi-sq = %.2f, df = %d, p = %s",
  kw_test$statistic, kw_test$parameter, format.pval(kw_test$p.value, digits = 3)
)

boxplot_quantile = ggplot(all_data, aes(x = quantile, y = sum_of_bits_combined_score, fill = quantile)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = kokako, guide = "none") +  
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = kw_label,
    hjust = -0.05, vjust = 2.0, size = 5.0
  ) +
  labs(
    x = "Quantile (ordered by degradation value)",
    y = "Sum of bits",
    title = "Sum of bits distribution across degradation quantiles (5' UTR)"
  ) +
  theme_bw(base_family = "arimo", base_size = 16)

# Display plots
print(scatter_plot)
print(boxplot_quantile)

# Save plots
ggsave("scatter.png", scatter_plot, width = 12, height = 4, dpi = 300)
ggsave("bitscore_by_quantile.png", boxplot_quantile, width = 12, height = 6, dpi = 300)

