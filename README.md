
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DrBoost

<!-- badges: start -->
<!-- badges: end -->

The goal of DrBoost is to provide **disjoint rule-based boosting** for
interpretable machine learning. It builds a piecewise-constant model as
a sum of short rules, where each observation is covered by exactly one
(non-overlapping) rule. This yields a globally interpretable
“one-rule-per-point” architecture, eliminating the ambiguity of
overlapping trees or rule sets.

## Installation

You can install the development version of DrBoost from
[GitHub](https://github.com/blind-contours/DrBoost) with:

``` r
# install.packages("pak")
pak::pak("blind-contours/DrBoost")

# Alternatively, if you use devtools you can do:
# devtools::install_github("blind-contours/DrBoost")
```

## Overview

DrBoost supports two main modes:

- **Logic-based for binary data**: discovers short AND/OR rules
  (disjoint).
- **Bounding-box for continuous features**: searches axis-aligned
  hyperrectangles using beam search.

Each chosen region is disjoint from the others, forming a partition of
the feature space. The coefficient on each rule (region) is solved by a
simple closed-form update (mean residual for continuous, or Newton step
for logistic). Our approach also offers early-stopping and optional
“exclude repeated failing region” logic to avoid adding harmful rules.

## Example 1: Binary Features and Logical Rules

Imagine you have a dataset with binary features X₁, X₂, …, Xₚ and a
binary outcome. Suppose the true label is governed by hidden logic like
(X₁=1 ∧ X₃=1) ∨ (X₂=0 ∧ X₅=1). We can try to recover that using DrBoost
in `logic=TRUE` mode.

``` r
library(DrBoost)

set.seed(123)
n <- 200

# Generate random binary features
X <- matrix(rbinom(n * 5, size=1, prob=0.3), nrow=n, ncol=5)
colnames(X) <- paste0("X", 1:5)

# True logic for Y: e.g. 1 if (X1=1 & X3=1) or (X2=0 & X5=1), else 0
true_logic <- function(row) {
  (row["X1"] == 1 && row["X3"] == 1) ||
    (row["X2"] == 0 && row["X5"] == 1)
}
Y <- apply(X, 1, true_logic)
Y <- as.integer(Y)  # Convert to 0/1

# Split data into train and validation
train_idx <- sample(seq_len(n), size = floor(0.7*n))
X_train <- X[train_idx, , drop=FALSE]
Y_train <- Y[train_idx]
X_val   <- X[-train_idx, , drop=FALSE]
Y_val   <- Y[-train_idx]

# Fit DrBoost in 'binary' mode with logic-based search
fit_logic <- disjoint_rule_boost(
  X_train, Y_train,
  X_val,   Y_val,
  outcome_type = "binary",
  logic       = TRUE,    # Use logic expansions
  max_iterations = 10,
  learning_rate  = 0.1,
  min_obs_pct    = 0.01,
  max_obs_frac   = 1.0,
  patience       = 2
)

# Inspect discovered rules
fit_logic$rule_list
# Each element has 'coverage', 'beta' (coefficient), and 'region_info$rule_string'
# indicating the logical expression.

# Evaluate on validation
pred_logodds_val <- fit_logic$f_val
pred_prob_val    <- 1 / (1 + exp(-pred_logodds_val))
cat("Final val loss =", fit_logic$val_loss_final, "\n")
```

DrBoost will display rules like:

    Iter 1 => coverage=..., val_loss=..., best_val=...
      rule: (X1==1 & X3==1)
    Iter 2 => coverage=..., val_loss=..., best_val=...
      rule: (X2==0 & X5==1)
    ...

Each observation ends up in exactly one rule region if they are
disjoint.

## Example 2: Continuous Features and Bounding Boxes

Now let’s show how DrBoost handles continuous features by searching for
axis-aligned boxes.

``` r
library(DrBoost)

set.seed(42)
n <- 300
X <- matrix(runif(n*2), nrow=n, ncol=2)
colnames(X) <- c("X1", "X2")

# Suppose a piecewise-constant true function:
f_true <- function(x1, x2) {
  if (x1 < 0.5 && x2 < 0.5) {
    return(2)
  } else if (x1 >= 0.5 && x2 < 0.5) {
    return(5)
  } else if (x1 < 0.5 && x2 >= 0.5) {
    return(3)
  } else {
    return(8)
  }
}

y_true <- apply(X, 1, function(r) f_true(r[1], r[2]))
y_obs  <- y_true + rnorm(n, sd=0.4)  # add some noise

# Split data
train_idx <- sample(seq_len(n), size=round(0.7*n))
X_train <- X[train_idx, , drop=FALSE]
Y_train <- y_obs[train_idx]
X_val   <- X[-train_idx, , drop=FALSE]
Y_val   <- y_obs[-train_idx]

# We need a 'thresholds_list' for bounding-box search:
build_thresholds_list <- function(X, num_thresholds=5) {
  p <- ncol(X)
  out <- matrix(NA_real_, nrow=num_thresholds, ncol=p)
  for (j in seq_len(p)) {
    out[, j] <- quantile(X[,j], probs=seq(0,1,length.out=num_thresholds))
  }
  out
}
thresholds_list <- build_thresholds_list(X_train, num_thresholds=5)

# Fit DrBoost
fit_box <- disjoint_rule_boost(
  X_train, Y_train,
  X_val,   Y_val,
  outcome_type     = "continuous",
  logic           = FALSE,           # bounding-box
  thresholds_list = thresholds_list,
  max_iterations  = 10,
  K1=50, K2=50, K3=50, K4=50,
  learning_rate = 0.1,
  patience = 3
)

# Inspect discovered bounding-box rules
fit_box$rule_list
# Each has a 'region_info$lowerBound' and 'region_info$upperBound',
# plus coverage and coefficient.

cat("Final validation MSE:", fit_box$val_loss_final, "\n")
```

You’ll see something like:

    Iter 1 => coverage=..., val_loss=..., best_val=...
      rule: (0 <= X1 <= 0.5) & (0 <= X2 <= 0.5)
    Iter 2 => coverage=..., val_loss=..., best_val=...
      rule: (0.5 < X1 <= 1) & (0 <= X2 <= 0.5)
    ...

indicating how DrBoost partitions the feature space with disjoint boxes
and assigns each its own coefficient.

## Additional Notes

- **Interpretability**: Because the rules do not overlap, every data
  point has exactly one “local explanation.”
- **Implementation**: The package uses a beam-search approach for
  bounding boxes, enumerating up to 4D subregions, and can optionally
  use OpenMP for parallel speedup.
- **Logic-based**: For binary data, setting `logic=TRUE` uses short
  AND/OR expansions (up to a specified literal depth).
- **Early Stopping**: The `patience` parameter halts training if we fail
  to improve validation loss for N consecutive steps.
- **Excluding failing regions**: An optional “exclude repeated failing
  region” mechanism avoids reselecting boxes that worsen validation.

We hope DrBoost helps you build interpretable, piecewise-constant
ensembles. For more details, please see the function documentation or
the project repository.
