# simulate_UCI_comparison.R
# ------------------------------------------------------------
# Example script to replicate the experiments described in the paper,
# using four UCI datasets. Compares Disjoint Rule Boosting (DRB)
# vs. CART, Random Forest, GBM, (optionally) RuleFit.
#
# Author: You
# ------------------------------------------------------------

library(DrBoost)         # Your package, containing disjoint_rule_boost()
library(rpart)           # For CART
library(randomForest)    # For Random Forest
library(gbm)             # For Gradient Boosting
# library(RuleFitPackage) # If you have a rulefit package available (optional)

library(MASS)            # For the Boston dataset (Housing)
# For Wine & Adult data, we might read from CSV or a known package
library(readr)

# 1) Helper function: train/val/test split
split_data <- function(X, Y, train_frac=0.70, val_frac=0.15, seed=123) {
  set.seed(seed)
  n <- nrow(X)
  idx <- sample(seq_len(n))
  
  train_cut  <- floor(train_frac * n)
  val_cut    <- floor((train_frac + val_frac) * n)
  
  train_idx <- idx[1:train_cut]
  val_idx   <- idx[(train_cut + 1):val_cut]
  test_idx  <- idx[(val_cut + 1):n]
  
  list(
    X_train = X[train_idx, , drop=FALSE],
    Y_train = Y[train_idx],
    X_val   = X[val_idx,   , drop=FALSE],
    Y_val   = Y[val_idx],
    X_test  = X[test_idx,  , drop=FALSE],
    Y_test  = Y[test_idx]
  )
}

# 2) Evaluate regression performance with RMSE
rmse <- function(pred, actual) {
  sqrt(mean((pred - actual)^2))
}

# 3) Evaluate classification performance with accuracy
accuracy <- function(pred_class, actual_class) {
  mean(pred_class == actual_class)
}

# 4) Run and compare for one dataset
#    data_type: "regression" or "classification"
#    returns a named vector or list with performance metrics
run_comparison <- function(X, Y, data_type=c("regression","classification"),
                           dataset_name="DatasetName") {
  data_type <- match.arg(data_type)
  
  # 70/15/15 split
  splitted <- split_data(X, Y, train_frac=0.70, val_frac=0.15, seed=1234)
  
  X_train <- splitted$X_train
  Y_train <- splitted$Y_train
  X_val   <- splitted$X_val
  Y_val   <- splitted$Y_val
  X_test  <- splitted$X_test
  Y_test  <- splitted$Y_test
  
  # 4a) Disjoint Rule Boosting
  # figure out feature types
  #   e.g. treat all numeric columns as 'continuous',
  #   unless they are 0/1 only => 'binary'
  # (optional) adjust for classification
  feat_types <- sapply(seq_len(ncol(X_train)), function(j) {
    vals <- unique(X_train[,j])
    if(all(vals %in% c(0,1))) "binary" else "continuous"
  })
  
  outcome_type <- ifelse(data_type=="classification","binary","continuous")
  
  # run DRB:
  drb_out <- disjoint_rule_boost(
    X_train=X_train, Y_train=Y_train,
    X_val=X_val,     Y_val=Y_val,
    feature_types=feat_types,
    outcome_type=outcome_type,
    max_iterations=50,   # or however many
    min_obs_pct=0.05,    # min coverage
    num_thresholds=20,   # number of threshold splits
    patience=3           # early stopping patience
  )
  
  # Evaluate on test
  #   For a final model, we can do:  predict(drb_out$final_model, newdata= <some DF> )
  if(data_type=="classification"){
    # create DF from X_test
    testDF <- as.data.frame(drb_out$rule_matrix_val)
    colnames(testDF) <- colnames(drb_out$rule_matrix_val)
    preds_link <- predict(drb_out$final_model, newdata=testDF, type="link")
    preds_prob <- 1 / (1 + exp(-preds_link))
    preds_class <- ifelse(preds_prob > 0.5, 1, 0)
    performance <- accuracy(preds_class, Y_test)
  } else {
    testDF <- as.data.frame(drb_out$rule_matrix_val)
    colnames(testDF) <- colnames(drb_out$rule_matrix_val)
    preds <- predict(drb_out$final_model, newdata=testDF)
    performance <- rmse(preds, Y_test)
  }
  
  # count how many final rules were used
  #   length(drb_out$regions) or the number of columns in rule_matrix_val
  n_rules_drb <- ncol(drb_out$rule_matrix_val)
  
  # 4b) CART
  # classification => rpart w/ method="class", else method="anova"
  if(data_type=="classification"){
    cart_method <- "class"
  } else {
    cart_method <- "anova"
  }
  cart_fit <- rpart(Y_train ~ ., data = data.frame(Y_train, X_train),
                    method=cart_method, control=rpart.control(cp=0.01))
  # prune or let rpart handle
  # predict on test
  if(data_type=="classification"){
    cart_probs <- predict(cart_fit, newdata=as.data.frame(X_test), type="prob")[,2]
    cart_preds <- ifelse(cart_probs>0.5,1,0)
    cart_perf <- accuracy(cart_preds, Y_test)
    # #Rules in a CART is typically #leaf-nodes?
    cart_rules <- length(unique(cart_fit$where))
  } else {
    cart_preds <- predict(cart_fit, newdata=as.data.frame(X_test))
    cart_perf <- rmse(cart_preds, Y_test)
    cart_rules <- length(unique(cart_fit$where))
  }
  
  # 4c) Random Forest
  rf_fit <- randomForest(x=X_train, y=Y_train, ntree=200)
  if(data_type=="classification"){
    rf_probs <- predict(rf_fit, newdata=X_test, type="prob")[,2]
    rf_preds <- ifelse(rf_probs>0.5,1,0)
    rf_perf <- accuracy(rf_preds, Y_test)
    rf_rules <- NA  # overlapping, so no well-defined single rule count
  } else {
    rf_preds <- predict(rf_fit, newdata=X_test)
    rf_perf <- rmse(rf_preds, Y_test)
    rf_rules <- NA
  }
  
  # 4d) GBM
  # classification => distribution="bernoulli", regression => distribution="gaussian"
  dist_gbm <- ifelse(data_type=="classification","bernoulli","gaussian")
  gbm_fit <- gbm::gbm.fit(
    x=X_train,
    y=Y_train,
    distribution=dist_gbm,
    n.trees=100,              # example
    interaction.depth=4,      # or whichever
    shrinkage=0.1,
    bag.fraction=1.0,         # no bagging here
    verbose=FALSE
  )
  best_iter <- gbm.perf(gbm_fit, method="OOB", plot.it=FALSE)
  if(data_type=="classification"){
    gbm_probs <- predict(gbm_fit, newdata=X_test, 
                         n.trees=best_iter, type="response")
    gbm_preds <- ifelse(gbm_probs>0.5,1,0)
    gbm_perf <- accuracy(gbm_preds, Y_test)
    gbm_rules <- NA
  } else {
    gbm_preds <- predict(gbm_fit, newdata=X_test, 
                         n.trees=best_iter, type="response")
    gbm_perf <- rmse(gbm_preds, Y_test)
    gbm_rules <- NA
  }
  
  # 4e) (Optional) RuleFit - if you have some function rulefit()
  # We'll skip real code as it depends on implementation.
  # Suppose something like:
  # rulefit_perf <- ...
  # rulefit_rules <- ...
  
  # Return results in a named vector or list
  out <- list(
    dataset = dataset_name,
    DRB_perf = performance,      # numeric (RMSE or ACC)
    DRB_rules= n_rules_drb,
    
    CART_perf= cart_perf,
    CART_rules= cart_rules,
    
    RF_perf  = rf_perf,
    RF_rules = rf_rules,
    
    GBM_perf = gbm_perf,
    GBM_rules= gbm_rules
    
    # If RuleFit included, add as well
    # RuleFit_perf = rulefit_perf,
    # RuleFit_rules= rulefit_rules
  )
  return(out)
}

# 5) Actually run for each dataset, store results in a data.frame
results_list <- list()

# A) Boston Housing (using MASS::Boston)
#     This is a regression problem, Y=medv
data("Boston", package="MASS")
X_boston <- as.matrix(Boston[, -which(colnames(Boston)=="medv")])
Y_boston <- Boston$medv

res_boston <- run_comparison(
  X=X_boston, Y=Y_boston, data_type="regression", 
  dataset_name="BostonHousing"
)
results_list[["BostonHousing"]] <- res_boston

# B) Wine Quality (red)
#   We can read from e.g. winequality-red.csv locally or from a known source
#   Suppose you have the CSV in "data/winequality-red.csv"
wine_red <- read_csv("data/winequality-red.csv")
# The official UCI data has columns like "fixed.acidity", "volatile.acidity", etc.
# The last column "quality" is the target
Y_wine <- wine_red$quality
X_wine <- as.matrix(wine_red[, setdiff(names(wine_red), "quality")])

res_wine <- run_comparison(
  X=X_wine, Y=Y_wine, data_type="regression",
  dataset_name="WineQualityRed"
)
results_list[["WineQualityRed"]] <- res_wine

# C) Adult dataset (classification)
#   Suppose you have "adult.data" loaded or use a package to get it
#   In typical usage, one must parse the CSV, convert >50K => 1, <=50K => 0, etc.
adult <- read_csv("data/adult.csv")  # Must match your local path
# You have to define X and Y properly, removing or encoding columns
# e.g. "income" is the target => convert to 0/1
adult$income_binary <- ifelse(adult$income == ">50K", 1, 0)
Y_adult <- adult$income_binary
# Suppose you exclude columns 'income' and something else
X_adult <- as.matrix(adult[, c("age","education.num","hours.per.week", "fnlwgt")]) # example

res_adult <- run_comparison(
  X=X_adult, Y=Y_adult, data_type="classification", 
  dataset_name="Adult"
)
results_list[["Adult"]] <- res_adult

# D) California Housing
#   There's a 'CaliforniaHousing' dataset from "mlbench" or "cali" from "AER" or we can read from CSV
#   Suppose we get from a CSV
cali <- read_csv("data/cal_housing.csv")
# The target might be 'median_house_value'
Y_cali <- cali$median_house_value
X_cali <- as.matrix(cali[, setdiff(names(cali), "median_house_value")])

res_cali <- run_comparison(
  X=X_cali, Y=Y_cali, data_type="regression",
  dataset_name="CaliforniaHousing"
)
results_list[["CaliforniaHousing"]] <- res_cali

# 6) Collect results into a data frame
df_results <- do.call(rbind, lapply(results_list, as.data.frame))

# 7) Print or save table
print(df_results)

# Optionally save to CSV
write.csv(df_results, "UCI_comparison_results.csv", row.names=FALSE)

message("Simulation script completed successfully!")
