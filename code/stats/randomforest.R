library(tidyverse)
library(ranger)
library(biomaRt)

# Load in FIMO data
fimo = read_tsv(c("fimo/stableq1fimo.tsv",
         "fimo/stableq2fimo.tsv",
         "fimo/stableq3fimo.tsv",
         "fimo/stableq4fimo.tsv",
         "fimo/stableq5fimo.tsv",
         "fimo/unstableq1fimo.tsv",
         "fimo/unstableq2fimo.tsv",
         "fimo/unstableq3fimo.tsv",
         "fimo/unstableq4fimo.tsv",
         "fimo/unstableq5fimo.tsv"),
         comment = "#")
# Load degradation spreadsheet
deg = read_csv("human_halflife_data_sorted.csv")

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

# Split data into train and test
set.seed(2)
train_index <- sample(1:nrow(model_data), 0.8 * nrow(model_data))
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Fit random forest 
# Handling step so R doesn't cry about the STREME motif names
predictor_cols <- setdiff(
  names(model_data),
  c("sequence_name", "ensembl_gene_id", "degradation_value")
)

set.seed(3)
rf_model <- ranger(
  data = train_data[, c(predictor_cols, "degradation_value")],
  dependent.variable.name = "degradation_value",
  importance = "permutation",
  num.trees = 1000
)

# Evaluate test set
predictions <- predict(rf_model, data = test_data[, predictor_cols])$predictions

rmse <- sqrt(mean((predictions - test_data$degradation_value)^2))
r_squared <- cor(predictions, test_data$degradation_value)^2
rmse
r_squared

# Extract and rank motif weights
motif_weights <- ranger::importance(rf_model) %>%
  enframe(name = "feature", value = "importance") %>%
  arrange(desc(importance))
View(motif_weights)
