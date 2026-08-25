# Genetic Algorithm–Based Feature Selection and Clustering Analysis

This repository provides a reproducible R workflow for genetic algorithm (GA)-based feature selection followed by hierarchical clustering and descriptive within-cohort evaluation.

Clinical disease-activity labels are used during the feature-selection stage through the GA fitness function. Therefore, the overall workflow is **label-guided**, although the final hierarchical clustering algorithm itself does not directly use the labels. The resulting ROC/AUC is an apparent within-cohort descriptive measure and should not be interpreted as an estimate of out-of-sample predictive performance.

## Repository contents

- `analysis.R` — runnable R script for GA feature selection, hierarchical clustering, centroid-distance score calculation, and descriptive ROC/AUC analysis.
- `synthetic_example_input.tsv` — synthetic samples × features matrix.
- `sample_annotation.tsv` — synthetic sample labels corresponding to the example feature matrix.
- `INPUT_DATA.md` — detailed input-data instructions.
- `sessionInfo.txt` — placeholder to be replaced by the `sessionInfo()` output from the R environment used for the manuscript analysis.
- `expected_output/` — representative output formats.

## Requirements

The manuscript analysis used R 4.4.1. The runnable script requires the following R packages:

```r
install.packages(c("GA", "pROC"))
```

Before public release, replace `sessionInfo.txt` with the actual output of `sessionInfo()` from the environment used to generate the manuscript results so that exact package versions are documented.

## Input data

See `INPUT_DATA.md` for the complete input specification. In brief, the feature matrix is a tab-delimited file with sample IDs in the first column and numerical features in all remaining columns. Clinical group labels are supplied in a separate `sample_annotation.tsv` file and matched by sample ID.

## Running the example

From the repository directory, run:

```bash
Rscript analysis.R
```

The script creates an `output/` directory automatically.

## Analysis workflow

1. Read the feature matrix and sample annotations.
2. Standardize each candidate feature subset.
3. Use hierarchical clustering with Ward's method and Euclidean distance.
4. Evaluate candidate subsets with a label-guided GA fitness function based on cluster-size and clinical-group-purity constraints.
5. Select the smallest feature subset among the optimal GA solutions.
6. Re-run hierarchical clustering using the selected features.
7. Calculate sample distances to the two cluster centroids and derive a continuous centroid-distance score.
8. Calculate a descriptive apparent within-cohort ROC/AUC.
9. Export selected features, cluster assignments, centroid-distance scores, the AUC summary, the GA optimization plot, and session information.

## Parameters

The draft script currently uses the values described in the revised manuscript Methods:

- minimum selected features: 2
- minimum cluster size: 12 samples
- minimum clinical-group purity per cluster: 0.70
- GA population size: 200
- maximum iterations: 200
- early stopping after 100 non-improving generations
- mutation probability: 0.2
- crossover probability: 0.5
- random seed: 42

These values must match the parameters actually used to generate the manuscript results before the repository is released.

## Output files

Running `analysis.R` generates:

- `output/selected_features.tsv`
- `output/cluster_assignments.tsv`
- `output/distance_scores.tsv`
- `output/auc_summary.tsv`
- `output/ga_optimization.pdf`
- `output/sessionInfo.txt`

Representative file formats are provided in `expected_output/`.

## Reproducibility note

The synthetic example data are provided only to demonstrate code execution and input/output structure. They do not reproduce the study's patient-level values or the manuscript's reported AUC.
