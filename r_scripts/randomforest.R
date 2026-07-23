library(tidyverse)
library(ranger)
library(biomaRt)
library(caret)
library(doParallel)
library(Manu)
library(sysfonts)
library(showtext)

# Font & colour setup 
font_add_google("Arimo", "arimo")
showtext_auto()
showtext_opts(dpi = 300)

kokako <- get_pal("Kokako")  

# Set region to calculate motif weights for
# Options: "all", "cds", "utr5", "utr3"
region <- "utr3"

# Load in FIMO data
fimo = read_tsv(c("plot_data/fimo/stable_scores/1stquan_scores/fimo.tsv",
                  "plot_data/fimo/stable_scores/2ndquan_scores/fimo.tsv",
                  "plot_data/fimo/stable_scores/3rdquan_scores/fimo.tsv",
                  "plot_data/fimo/stable_scores/4thquan_scores/fimo.tsv",
                  "plot_data/fimo/stable_scores/5thquan_scores/fimo.tsv",
                  "plot_data/fimo/unstable_scores/1stquan_scores/fimo.tsv",
                  "plot_data/fimo/unstable_scores/2ndquan_scores/fimo.tsv",
                  "plot_data/fimo/unstable_scores/3rdquan_scores/fimo.tsv",
                  "plot_data/fimo/unstable_scores/4thquan_scores/fimo.tsv",
                  "plot_data/fimo/unstable_scores/5thquan_scores/fimo.tsv"),
         comment = "#")

# Load degradation spreadsheet
deg = read_csv("plot_data/human_halflife_data_sorted.csv")

# Load transcript region boundaries CSV
boundaries = read_csv("plot_data/boundaries.csv")

# Filter FIMO hits to select those in chosen region
if (region != "all") {
  joined <- fimo %>%
    inner_join(boundaries, by = c("sequence_name" = "sequence_id")) %>%
    mutate(hit_mid = (start + stop) / 2)
  
  fimo <- switch(region,
    "cds" = joined %>% filter(hit_mid > cds_start & hit_mid <= cds_end),
    "utr5" = joined %>% filter(hit_mid > utr5_start & hit_mid <= utr5_end),
    "utr3" = joined %>% filter(hit_mid > utr3_start & hit_mid <= utr3_end)
    ) %>%
    dplyr::select(-hit_mid, -seq_len, -has_orf, -utr5_start, -utr5_end, -utr5_len,
                  -cds_start, -cds_end, -cds_len, -utr3_start, -utr3_end, -utr3_len)
}

# Collapse FIMO hits into transcript-motif pair scores
motif_scores <- fimo %>%
  mutate(neglog10_q = -log10(`q-value`)) %>%
  group_by(sequence_name, motif_id) %>%
  summarise(motif_score = sum(neglog10_q), .groups = "drop")

# Reshape into matrix
feature_matrix <- motif_scores %>%
  pivot_wider(names_from = motif_id, values_from = motif_score, values_fill = 0)

# Convert transcript IDs to gene IDs
# Strip verison numbers from transcript IDs 
feature_matrix <- feature_matrix %>% mutate(sequence_name = sub("\\..*", "", sequence_name))

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

id_map <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
  filters = "ensembl_transcript_id",
  values = feature_matrix$sequence_name,
  mart = ensembl
)

feature_matrix <- feature_matrix %>%
  left_join(id_map, by = c("sequence_name" = "ensembl_transcript_id"))

# Attach degradation values
model_data <- feature_matrix %>%
  inner_join(deg %>% dplyr::select(gene_id, degradation_value),
             by = c("ensembl_gene_id" = "gene_id"))

# Add random number column
set.seed(1)
model_data$random_numbers <- runif(nrow(model_data))

# Identify predictor columns
predictor_cols <- setdiff(
  names(model_data),
  c("sequence_name", "ensembl_gene_id", "degradation_value")
)

# Split data into train and test
set.seed(2)
train_index <- createDataPartition(model_data$degradation_value, p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Fit random forest 
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 1,
  savePredictions = "final",
)

p <- length(predictor_cols)
rf_grid <- expand.grid(
  mtry = 20,
  splitrule = "variance",
  min.node.size = 1
)

cl <- makePSOCKcluster(parallel::detectCores() - 1)   
registerDoParallel(cl)

set.seed(3)
rf_fit <- train(
  x = train_data[, predictor_cols],
  y = train_data$degradation_value,
  method = "ranger",
  importance = "permutation",
  num.trees = 500,
  trControl = ctrl,
  tuneGrid = rf_grid,
  metric = "RMSE",
  num.threads = 1
)

stopCluster(cl)         
registerDoSEQ()

rf_fit
# Only need to run if you're tuning model parameters
plot(rf_fit)

# Evaluate test set
predictions <- predict(rf_fit, newdata = test_data[, predictor_cols])

test_performance <- postResample(pred = predictions, obs = test_data$degradation_value)
test_performance

# Cross-validate performance from training
cv_performance <- rf_fit$resample
cat("CV RMSE: ", round(mean(cv_performance$RMSE), 4),
    " +/- ", round(sd(cv_performance$RMSE), 4), "\n")
cat("CV Rsquared: ", round(mean(cv_performance$Rsquared), 4),
    " +/- ", round(sd(cv_performance$Rsquared), 4), "\n")

# Extract and rank motif weights
motif_weights <- varImp(rf_fit, scale = FALSE)$importance %>%
  rownames_to_column("feature") %>%
  rename(importance = Overall) %>%
  arrange(desc(importance))

random_benchmark <- motif_weights %>%
  filter(feature == "random_numbers") %>%
  pull(importance)

motif_weights <- motif_weights %>%
  mutate(above_noise_floor = importance > random_benchmark) %>%
  filter(feature != "random_numbers")

# Rename motifs
parse_streme_ids <- function(path, sign) {
  motif_lines <- grep("^MOTIF", readLines(path), value = TRUE)
  ids <- sub("^MOTIF\\s+(\\S+).*", "\\1", motif_lines)
  tibble(
    motif_id = ids,
    motif_stripped = sub("^[0-9]+-", "", ids),
    sign = sign
  )
}

motif_key <- bind_rows(
  parse_streme_ids("plot_data/streme_stable.txt", "+"),
  parse_streme_ids("plot_data/streme_unstable.txt", "-")
)

motif_weights <- motif_weights %>%
  left_join(motif_key, by = c("feature" = "motif_id")) %>%
  mutate(motif_label = if_else(is.na(sign), feature, paste0(motif_stripped, sign))) %>%
  dplyr::select(-motif_stripped, -sign)

print(motif_weights)

# Make output directory
output_dir <- file.path("rf_results", region)
dir.create(output_dir, recursive = TRUE)

# Save random forest model
saveRDS(rf_fit, file.path(output_dir, "rf_fit.rds"))

# Save results
write_csv(as.data.frame(t(test_performance)), file.path(output_dir, "test_performance.csv"))
write_csv(motif_weights, file.path(output_dir, "motif_weights.csv"))
cv_summary <- tibble(
  Resample = c("Mean", "SD"),
  RMSE = c(mean(cv_performance$RMSE), sd(cv_performance$RMSE)),
  Rsquared = c(mean(cv_performance$Rsquared), sd(cv_performance$Rsquared)),
  MAE = c(mean(cv_performance$MAE), sd(cv_performance$MAE))
)
cv_performance_out <- bind_rows(cv_performance, cv_summary)
write_csv(cv_performance_out, file.path(output_dir, "cv_performance.csv"))

## Plotting functions
# Obs vs prd plot
obs_prd_plot <- ggplot(data.frame(observed = test_data$degradation_value, predicted = predictions),
                       aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5, colour = kokako[1]) +
  geom_abline(intercept = 0, slope = 1, color = "grey40", linetype = "dashed") +
  labs(title = "Observed vs Predicted Degradation Values",
       x = "Predicted", y = "Observed") +
  theme_minimal(base_family = "arimo", base_size = 20)

obs_prd_plot
ggsave(file.path(output_dir, "obs_vs_prd_plot.png"), obs_prd_plot, width = 6, height = 5)

# Residual plot
residual_plot <- ggplot(data.frame(predicted = predictions,
                                   residual = test_data$degradation_value - predictions),
                        aes(x = predicted, y = residual)) +
  geom_point(alpha = 0.5, colour = kokako[1]) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
  labs(title = "Residuals vs Predicted Values",
       x = "Predicted", y = "Residual") +
  theme_minimal(base_family = "arimo", base_size = 20)

residual_plot
ggsave(file.path(output_dir, "residuals_plot.png"), residual_plot, width = 6, height = 5)

# Variable importance bar chart
var_imp_plot <- ggplot(motif_weights, aes(x = reorder(motif_label, importance), y = importance)) +
  geom_col(fill = kokako[3]) +
  coord_flip() +
  labs(title = "Motif Importance", x = "Motif occurrences across 3' UTR", y = "Importance") +
  theme_minimal(base_family = "arimo", base_size = 16)

var_imp_plot
ggsave(file.path(output_dir, "motif_importance_plot.png"), var_imp_plot, width = 7, height = 8)

# CV fold variability plot
cv_long <- cv_performance %>%
  dplyr::select(Resample, RMSE, Rsquared, MAE) %>%
  pivot_longer(cols = c(RMSE, Rsquared, MAE), names_to = "Metric", values_to = "Value")

cv_plot <- ggplot(cv_long, aes(x = Metric, y = Value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, colour = kokako[2]) +
  facet_wrap(~ Metric, scales = "free") +
  labs(title = "Cross-Validation Performance Across Folds", x = NULL, y = NULL) +
  theme_minimal(base_family = "arimo", base_size = 20)

cv_plot
ggsave(file.path(output_dir, "cv_fold_variability_plot.png"), cv_plot, width = 8, height = 4)

