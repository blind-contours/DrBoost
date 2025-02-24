# comparison_simulation.R
# ----------------------------------------------------
# Runs Disjoint Rule Boosting vs. CART, RF, GBM, XGBoost
# on five small UCI regression datasets (like Boston housing, etc.).

# We'll assume your new DRB code + predict_disjoint_rule_boost() are loaded.

DRBOOST_PARAMS <- list(
  max_iterations = 20,
  num_thresholds = 18,
  min_obs_pct    = 0.01,
  max_obs_frac = 0.25,
  patience       = 10,
  K1 = 300,
  K2 = 300,
  K3 = 300,
  K4 = 300,
  learning_rate = 9.0
)

library(readr)
library(rpart)
library(randomForest)
library(gbm)
library(xgboost)
library(devtools)
library(here)

set.seed(42)

mse <- function(pred, actual) {
  mean((pred - actual)^2)
}

split_data <- function(X, Y, seed=123, train_frac=0.70, val_frac=0.15) {
  set.seed(seed)
  n <- nrow(X)
  idx <- sample(seq_len(n))
  train_cut <- floor(train_frac * n)
  val_cut   <- floor((train_frac + val_frac) * n)
  
  train_idx <- idx[1:train_cut]
  val_idx   <- idx[(train_cut+1):val_cut]
  test_idx  <- idx[(val_cut+1):n]
  
  list(
    X_train = X[train_idx, , drop=FALSE],
    Y_train = Y[train_idx],
    X_val   = X[val_idx,   , drop=FALSE],
    Y_val   = Y[val_idx],
    X_test  = X[test_idx,  , drop=FALSE],
    Y_test  = Y[test_idx]
  )
}


run_comparison_regression <- function(dataset_name, X, Y,
                                      drboost_params=DRBOOST_PARAMS) {
  cat("\n------------------------------\n")
  cat("DATASET:", dataset_name, "\n")
  cat("------------------------------\n")
  
  # 1) Split data
  splitted <- split_data(X, Y, seed=2024)
  X_train <- splitted$X_train; Y_train <- splitted$Y_train
  X_val   <- splitted$X_val;   Y_val   <- splitted$Y_val
  X_test  <- splitted$X_test;  Y_test  <- splitted$Y_test
  
  # rename columns
  colnames(X_train) <- sprintf("X%d", seq_len(ncol(X_train)))
  colnames(X_val)   <- sprintf("X%d", seq_len(ncol(X_val)))
  colnames(X_test)  <- sprintf("X%d", seq_len(ncol(X_test)))
  
  #------------------------------
  # 2) Disjoint Rule Boosting
  #------------------------------
  thresholds_list <- build_thresholds_list(
    X_train,
    num_thresholds = drboost_params$num_thresholds
  )
  
  # *** Using your new DRB approach that returns (init, rule_list, etc.)
  drb_out <- disjoint_rule_boost(
    X_train       = X_train,
    Y_train       = Y_train,
    X_val         = X_val,
    Y_val         = Y_val,
    outcome_type  = "continuous",  # regression
    max_iterations= drboost_params$max_iterations,
    min_obs_pct   = drboost_params$min_obs_pct,
    max_obs_frac = drboost_params$max_obs_frac,
    thresholds_list= thresholds_list,
    patience      = drboost_params$patience,
    K1            = drboost_params$K1,
    K2            = drboost_params$K2,
    K3            = drboost_params$K3,
    K4            = drboost_params$K4,
    # possibly pass learning_rate=..., second_order_logistic=FALSE, etc.
    featureNames  = colnames(X_train),
    learning_rate = drboost_params$learning_rate
  )
  
  # Now predict on test data
  # because we have a "continuous" outcome, we just do "type='link'" or no arg
  test_pred_drb <- predict_disjoint_rule_boost(
    object       = drb_out,
    X_new        = X_test,
    outcome_type = "continuous", 
    logic        = FALSE,         # bounding-box
    type         = "link"         # raw numeric
  )
  drb_mse <- mse(test_pred_drb, Y_test)
  
  # number of rules discovered
  drb_rules_ct <- length(drb_out$rule_list)
  
  if (drb_rules_ct > 0) {
    cat("** First discovered rule for DRB:\n")
    cat("   Rule:", drb_out$rule_list[[1]]$region_info$rule_string, "\n")
  }
  
  #------------------------------
  # 3) CART
  #------------------------------
  cart_fit <- rpart(
    Y_train ~ .,
    data = data.frame(Y_train, X_train),
    method="anova",
    control = rpart.control(cp=0.01, maxdepth = 4)
  )
  cart_pred_test <- predict(cart_fit, newdata = data.frame(X_test))
  cart_mse <- mse(cart_pred_test, Y_test)
  cart_leaves <- length(unique(cart_fit$where))
  
  #------------------------------
  # 4) Random Forest
  #------------------------------
  rf_fit <- randomForest(x=X_train, y=Y_train, ntree=200)
  rf_pred_test <- predict(rf_fit, newdata=X_test)
  rf_mse <- mse(rf_pred_test, Y_test)
  
  #------------------------------
  # 5) GBM
  #------------------------------
  gbm_fit <- gbm::gbm.fit(
    x = X_train, y = Y_train,
    distribution="gaussian",
    n.trees=100, interaction.depth=4, shrinkage=0.1,
    bag.fraction=0.5, verbose=FALSE
  )
  best_iter <- gbm.perf(gbm_fit, method="OOB", plot.it=FALSE)
  gbm_pred_test <- predict(gbm_fit, newdata=X_test, n.trees=best_iter)
  gbm_mse <- mse(gbm_pred_test, Y_test)
  
  #------------------------------
  # 6) XGBoost
  #------------------------------
  dtrain <- xgb.DMatrix(data=X_train, label=Y_train)
  dtest  <- xgb.DMatrix(data=X_test,  label=Y_test)
  xgb_params <- list(objective="reg:squarederror", max_depth=4, eta=0.1, subsample=0.5)
  xgb_model <- xgb.train(
    params   = xgb_params,
    data     = dtrain,
    nrounds  = 200,
    watchlist= list(train=dtrain),
    verbose  = 0
  )
  xgb_pred_test <- predict(xgb_model, newdata=X_test)
  xgb_mse <- mse(xgb_pred_test, Y_test)
  
  #------------------------------
  # 7) Return summary
  #------------------------------
  out <- list(
    dataset      = dataset_name,
    DRB_MSE      = drb_mse,
    DRB_rules    = drb_rules_ct,
    CART_MSE     = cart_mse,
    CART_leaves  = cart_leaves,
    RF_MSE       = rf_mse,
    GBM_MSE      = gbm_mse,
    XGB_MSE      = xgb_mse,
    
    # keep the DRB object, or we can store just init & rule_list
    DRB_out      = drb_out,
    XGB_model    = xgb_model
  )
  return(out)
}


#----------------------------------------
# Main Execution
#----------------------------------------
all_results <- list()
data_files <- c("boston.csv", "concrete.csv", "energy.csv", "wine.csv", "yacht.csv")

for(fname in data_files) {
  fpath <- here(file.path("sandbox/data", fname))
  df <- read_csv(fpath, col_types = cols(.default=col_double()))
  
  # last column is outcome
  Y <- df[[ncol(df)]]
  X <- df[, -ncol(df), drop=FALSE]
  
  X_mat <- as.matrix(X)
  dataset_label <- sub(".csv$","", fname)
  
  res <- run_comparison_regression(
    dataset_name = dataset_label,
    X = X_mat,
    Y = Y,
    drboost_params = DRBOOST_PARAMS
  )
  
  all_results[[dataset_label]] <- res
}

# Combine summary
df_results <- do.call(rbind, lapply(all_results, function(x) {
  as.data.frame(x[setdiff(names(x),
                          c("DRB_out","XGB_model"))])
}))
print(df_results)
cat("\nAll done. MSE etc. in df_results.\n")
