# Genetic Algorithm–Based Feature Selection and Clustering Analysis
#
# This workflow is label-guided at the feature-selection stage because clinical
# disease-activity labels enter the GA fitness function. The final hierarchical
# clustering algorithm itself is unsupervised.

suppressPackageStartupMessages({
  library(GA)
  library(pROC)
})

set.seed(42)

# -----------------------------
# Parameters
# -----------------------------
min_selected_features <- 2
min_cluster_size <- 12
purity_threshold <- 0.70

# -----------------------------
# Input
# -----------------------------
data <- read.table(
  "synthetic_example_input.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

annotation <- read.table(
  "sample_annotation.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (anyDuplicated(rownames(data))) {
  stop("Duplicated sample IDs were detected in the feature matrix.")
}
if (anyDuplicated(annotation$sample_id)) {
  stop("Duplicated sample IDs were detected in sample_annotation.tsv.")
}
if (!all(rownames(data) %in% annotation$sample_id)) {
  stop("Some samples in the feature matrix do not have annotations.")
}

annotation <- annotation[match(rownames(data), annotation$sample_id), , drop = FALSE]
annotation_row <- factor(annotation$status, levels = c("stable", "active"))

if (any(is.na(annotation_row))) {
  stop("Status must be either 'stable' or 'active' for every sample.")
}

is_numeric <- vapply(data, is.numeric, logical(1))
if (!all(is_numeric)) {
  stop(
    "All feature columns must be numeric. Non-numeric columns: ",
    paste(colnames(data)[!is_numeric], collapse = ", ")
  )
}
if (anyNA(data)) {
  stop("Missing feature values were detected. Please resolve them before analysis.")
}

feature_names <- colnames(data)
num_features <- ncol(data)

# -----------------------------
# GA fitness function
# -----------------------------
evaluate_subset <- function(chromosome) {
  selected <- which(chromosome == 1)

  if (length(selected) < min_selected_features) {
    return(0)
  }

  signature_values <- scale(data[, selected, drop = FALSE])

  if (any(!is.finite(signature_values))) {
    return(0)
  }

  dend <- hclust(dist(signature_values), method = "ward.D")
  cutree_2 <- cutree(dend, k = 2)

  # Require sufficiently large clusters.
  if (min(table(cutree_2)) < min_cluster_size) {
    return(0)
  }

  tab <- table(status = annotation_row, cluster = cutree_2)

  # Clinical-group purity within each cluster.
  cluster_purity <- apply(tab, 2, function(x) max(x) / sum(x))

  if (!all(cluster_purity >= purity_threshold)) {
    return(0)
  }

  # Cluster numbers are arbitrary, so evaluate both possible mappings.
  score_mapping_1 <- tab[1, 1] + tab[2, 2]
  score_mapping_2 <- tab[1, 2] + tab[2, 1]

  max(score_mapping_1, score_mapping_2)
}

# -----------------------------
# Run GA
# -----------------------------
ga_result <- ga(
  type = "binary",
  fitness = evaluate_subset,
  nBits = num_features,
  popSize = 200,
  maxiter = 200,
  run = 100,
  pmutation = 0.2,
  pcrossover = 0.5,
  seed = 42
)

print(summary(ga_result))

if (!dir.exists("output")) {
  dir.create("output")
}

pdf("output/ga_optimization.pdf")
plot(ga_result)
dev.off()

# -----------------------------
# Select the smallest optimal solution
# -----------------------------
best_solution <- ga_result@solution
if (is.null(dim(best_solution))) {
  best_solution <- matrix(best_solution, nrow = 1)
}

solution_sizes <- rowSums(best_solution)
best_index <- which.min(solution_sizes)
chosen_solution <- best_solution[best_index, ]
selected_features <- feature_names[chosen_solution == 1]

write.table(
  data.frame(feature = selected_features),
  file = "output/selected_features.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Final hierarchical clustering
# -----------------------------
selected_data <- scale(data[, selected_features, drop = FALSE])
final_dend <- hclust(dist(selected_data), method = "ward.D")
final_cluster <- cutree(final_dend, k = 2)

cluster_table <- table(status = annotation_row, cluster = final_cluster)

# Orient cluster numbering to clinical groups only for descriptive interpretation.
score_alignment_1 <- cluster_table["stable", "1"] + cluster_table["active", "2"]
score_alignment_2 <- cluster_table["stable", "2"] + cluster_table["active", "1"]

if (score_alignment_1 >= score_alignment_2) {
  stable_cluster <- 1
  active_cluster <- 2
} else {
  stable_cluster <- 2
  active_cluster <- 1
}

cluster_assignments <- data.frame(
  sample_id = rownames(data),
  status = as.character(annotation_row),
  cluster = final_cluster
)

write.table(
  cluster_assignments,
  file = "output/cluster_assignments.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Cluster-centroid distance score
# -----------------------------
active_centroid <- colMeans(
  selected_data[final_cluster == active_cluster, , drop = FALSE]
)
stable_centroid <- colMeans(
  selected_data[final_cluster == stable_cluster, , drop = FALSE]
)

euclidean_distance <- function(x, centroid) {
  sqrt(sum((x - centroid)^2))
}

distance_to_active <- apply(
  selected_data,
  1,
  euclidean_distance,
  centroid = active_centroid
)

distance_to_stable <- apply(
  selected_data,
  1,
  euclidean_distance,
  centroid = stable_centroid
)

# Larger values indicate greater proximity to the active-cluster centroid.
separation_score <- distance_to_stable - distance_to_active

distance_scores <- data.frame(
  sample_id = rownames(data),
  status = as.character(annotation_row),
  cluster = final_cluster,
  distance_to_active_centroid = distance_to_active,
  distance_to_stable_centroid = distance_to_stable,
  separation_score = separation_score
)

write.table(
  distance_scores,
  file = "output/distance_scores.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Descriptive apparent within-cohort ROC/AUC
# -----------------------------
roc_result <- roc(
  response = annotation_row,
  predictor = separation_score,
  levels = c("stable", "active"),
  direction = "<",
  quiet = TRUE
)

auc_value <- as.numeric(auc(roc_result))

write.table(
  data.frame(
    metric = "apparent_within_cohort_AUC",
    value = auc_value
  ),
  file = "output/auc_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Selected features:\n")
print(selected_features)
cat("\nApparent within-cohort AUC:", auc_value, "\n")

# Record R and package versions used for this run.
capture.output(sessionInfo(), file = "output/sessionInfo.txt")
