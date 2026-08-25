# Input data instructions

The analysis uses two tab-delimited input files: a feature matrix and a sample-annotation file.

## 1. Feature matrix: `synthetic_example_input.tsv`

Rows represent individual samples and columns represent numerical features used for genetic algorithm (GA)-based feature selection.

- The first column contains unique sample identifiers.
- The first field of the header row is left empty.
- All remaining columns are numerical features.
- Feature names must be unique.
- Missing values are not allowed in the example workflow.

Example format:

```text
\tfeature1\tfeature2\tfeature3
donor1\t11.2\t0.0920392281057903\t5.35
donor2\t45.42\t0.360581485898529\t11.53
donor3\t24.02\t0.12642462376519\t12.77
```

The file is read in R as:

```r
data <- read.table(
  "synthetic_example_input.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
```

## 2. Sample annotation: `sample_annotation.tsv`

The annotation file contains two columns:

- `sample_id`: unique sample identifier matching the row names of the feature matrix.
- `status`: clinical disease-activity label. In this example, the allowed values are `active` and `stable`.

Example:

```text
sample_id\tstatus
donor1\tactive
donor2\tactive
donor18\tstable
```

The script matches annotations to the feature matrix by `sample_id`, so the files do not need to be manually reordered before analysis. Duplicate sample IDs, missing annotations, non-numeric feature columns, and missing feature values will trigger an error.

## Synthetic data note

The provided `synthetic_example_input.tsv` and `sample_annotation.tsv` are synthetic and are supplied only to demonstrate the required file structure and execution of the analysis workflow. They do not contain patient-level study data and are not intended to reproduce the numerical results reported in the manuscript.
