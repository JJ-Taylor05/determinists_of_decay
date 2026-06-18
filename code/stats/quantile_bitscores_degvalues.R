library(tidyverse)
library(boot)

# Read in quantile files
quantile_files = c("q1" = "merged_combined_scores_q1.tsv",
                   "q2" = "merged_combined_scores_q2.tsv",
                   "q3" = "merged_combined_scores_q3.tsv",
                   "q4" = "merged_combined_scores_q4.tsv",
                   "q5" = "merged_combined_scores_q5.tsv")

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
  geom_smooth(method = "loess", se = FALSE, colour = "plum") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = sprintf("rho = %.2f\n[%.2f, %.2f]\nn = %d", 
                    overall_result$rho, 
                    overall_result$ci_lower, 
                    overall_result$ci_upper,
                    overall_result$n),
    hjust = -0.1, vjust = 1.2, size = 3.5
  ) +
  labs(
    x = "Degradation value",
    y = "Sum of bits",
    colour = "Quantile",
    title = "Degradation values vs sum of bits score"
  ) +
  theme_bw()

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

boxplot_quantile = ggplot(all_data, aes(x = quantile, y = sum_of_bits_combined_score)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "lightskyblue") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = kw_label,
    hjust = -0.05, vjust = 2.0, size = 3.5
  ) +
  labs(
    x = "Quantile (ordered by degradation value)",
    y = "Sum of bits",
    title = "Sum of bits distribution across degradation quantiles"
  ) +
  theme_bw()

# Display plots
print(scatter_plot)
print(boxplot_quantile)
