# comparison_simulation.R
# ----------------------------------------------------
# Runs Disjoint Rule Boosting vs CART, RF, GBM, XGBoost
# on five small UCI regression datasets (Boston, Concrete, Energy, Wine, Yacht),
# using local_model_type="linear" for piecewise-linear expansions
# in the DRB model, with leftover region fitting.

# 1) Adjusted DRBOOST_PARAMS for thorough coverage & leftover refinement
DRBOOST_PARAMS <- list(
  # pass 1 config
  max_pass1_iter = 100,  # more iterations => more boxes if beneficial
  min_obs_pct    = 0.05, # allow smaller coverage
  max_obs_frac   = 1.0,   # up to 100% coverage if it helps
  patience       = 10,    # more relaxed => picks more boxes if incremental improvement
  
  # expansions & learning rate
  num_thresholds = 10,    # how many threshold splits per feature
  K1 = 500,
  K2 = 4000,
  K3 = 300,
  K4 = 200,
  learning_rate  = 1.0,
  
  # coverage overlap
  max_overlap    = 1,     # 1 => strictly disjoint (here we allow overlap up to 2)
  
  # leftover region refinement
  refine_leftover= TRUE
)

library(readr)
library(rpart)
library(randomForest)
library(gbm)
library(xgboost)
library(devtools)
library(here)

set.seed(42)

# Basic MSE function
mse <- function(pred, actual) {
  mean((pred - actual)^2)
}

#------------------------------------------------------------------
# Utility: split data into train/val/test
#------------------------------------------------------------------
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

#------------------------------------------------------------------
# Build thresholds for DRB
# Typically done via quantiles. We'll store them in a matrix
# where row=threshold, col=feature.
#------------------------------------------------------------------
build_thresholds_list <- function(X, num_thresholds=18) {
  p <- ncol(X)
  out <- matrix(NA_real_, nrow=num_thresholds, ncol=p)
  for(j in seq_len(p)) {
    xj <- X[,j]
    qs <- quantile(xj, probs=seq(0,1,length.out=num_thresholds), na.rm=TRUE)
    out[,j] <- qs
  }
  out
}

#------------------------------------------------------------------
# Main function to run experiments for one dataset
#------------------------------------------------------------------
run_comparison_regression <- function(dataset_name, X, Y,
                                      drboost_params=DRBOOST_PARAMS) {
  cat("\n------------------------------\n")
  cat("DATASET:", dataset_name, "\n")
  cat("------------------------------\n")
  
  # 1) Split data
  splitted <- split_data(X, Y, seed=2024, train_frac=0.70, val_frac=0.15)
  X_train <- splitted$X_train; Y_train <- splitted$Y_train
  X_val   <- splitted$X_val;   Y_val   <- splitted$Y_val
  X_test  <- splitted$X_test;  Y_test  <- splitted$Y_test
  
  # Force rename columns => "X1","X2",..., for consistency with some ML packages
  colnames(X_train) <- paste0("X", seq_len(ncol(X_train)))
  colnames(X_val)   <- paste0("X", seq_len(ncol(X_val)))
  colnames(X_test)  <- paste0("X", seq_len(ncol(X_test)))
  
  # 2) Build thresholds for DRB
  thresholds_list <- build_thresholds_list(
    X_train,
    num_thresholds = drboost_params$num_thresholds
  )
  
  # 3) Disjoint Rule Boosting (piecewise-linear expansions)
  cat("Running Disjoint Rule Boosting...\n")
  drb_out <- disjoint_rule_boost(
    X_train       = X_train,
    Y_train       = Y_train,
    X_val         = X_val,
    Y_val         = Y_val,
    outcome_type  = "continuous",
    
    # main parameters from drboost_params
    min_obs_pct   = drboost_params$min_obs_pct,
    max_obs_frac  = drboost_params$max_obs_frac,
    thresholds_list = thresholds_list,
    patience      = drboost_params$patience,
    K1            = drboost_params$K1,
    K2            = drboost_params$K2,
    K3            = drboost_params$K3,
    K4            = drboost_params$K4,
    learning_rate = drboost_params$learning_rate,
    
    # pass 1 iteration limit
    max_pass1_iter = drboost_params$max_pass1_iter,
    
    # refine dimension
    max_dim_refine = 4,
    
    # coverage overlap
    max_overlap    = drboost_params$max_overlap,
    
    # leftover region refinement
    refine_leftover= drboost_params$refine_leftover
  )
  
  cat("Predicting with DRB...\n")
  test_pred_drb <- predict_disjoint_rule_boost(
    object       = drb_out,
    X_new        = X_test,
    outcome_type = "continuous",
    type         = "link"  # raw numeric predictions
  )
  drb_mse <- mse(test_pred_drb, Y_test)
  drb_rules_ct <- length(drb_out$rule_list)
  
  cat(sprintf("DRB => MSE=%.4f, rules=%d\n", drb_mse, drb_rules_ct))
  
  # Optionally show summary of region rules
  cat("\n*** Disjoint Rule Boosting Summary of Region Rules ***\n")
  # If you have a summary.disjoint_rule_boost() method:
  drb_summary <- summary(drb_out, X_train=X_train, Y_train=Y_train)
  
  #------------------------------
  # 4) CART
  #------------------------------
  cat("\nFitting CART...\n")
  cart_fit <- rpart(
    Y_train ~ .,
    data = data.frame(Y_train, X_train),
    method="anova",
    control = rpart.control(cp=0.01, maxdepth=4)
  )
  cart_pred_test <- predict(cart_fit, newdata=data.frame(X_test))
  cart_mse <- mse(cart_pred_test, Y_test)
  cart_leaves <- length(unique(cart_fit$where))
  cat(sprintf("CART => MSE=%.4f, leaves=%d\n", cart_mse, cart_leaves))
  
  #------------------------------
  # 5) Random Forest
  #------------------------------
  cat("Fitting Random Forest...\n")
  rf_fit <- randomForest(x=X_train, y=Y_train, ntree=200)
  rf_pred_test <- predict(rf_fit, newdata=X_test)
  rf_mse <- mse(rf_pred_test, Y_test)
  cat(sprintf("RF => MSE=%.4f\n", rf_mse))
  
  #------------------------------
  # 6) GBM
  #------------------------------
  cat("Fitting GBM...\n")
  gbm_fit <- gbm::gbm.fit(
    x = as.data.frame(X_train), y = Y_train,
    distribution="gaussian",
    n.trees=100, interaction.depth=4, shrinkage=0.1,
    bag.fraction=0.5, verbose=FALSE
  )
  best_iter <- gbm.perf(gbm_fit, method="OOB", plot.it=FALSE)
  gbm_pred_test <- predict(gbm_fit, newdata=X_test, n.trees=best_iter)
  gbm_mse <- mse(gbm_pred_test, Y_test)
  cat(sprintf("GBM => MSE=%.4f\n", gbm_mse))
  
  #------------------------------
  # 7) XGBoost
  #------------------------------
  cat("Fitting XGBoost...\n")
  dtrain <- xgb.DMatrix(data=as.matrix(X_train), label=Y_train)
  xgb_params <- list(objective="reg:squarederror", max_depth=4, eta=0.1, subsample=0.5)
  xgb_model <- xgb.train(
    params   = xgb_params,
    data     = dtrain,
    nrounds  = 200,
    verbose  = 0
  )
  xgb_pred_test <- predict(xgb_model, newdata=as.matrix(X_test))
  xgb_mse <- mse(xgb_pred_test, Y_test)
  cat(sprintf("XGB => MSE=%.4f\n", xgb_mse))
  
  #------------------------------
  # 8) Summarize & Return
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
    
    DRB_out      = drb_out,
    DRB_summary  = drb_summary  # store the summary if needed
  )
  return(out)
}

#----------------------------------------
# Main Execution over multiple data files
#----------------------------------------
all_results <- list()
data_files <- c("boston.csv", "concrete.csv", "energy.csv", "wine.csv", "yacht.csv")

for(fname in data_files) {
  fpath <- here::here(file.path("sandbox/data", fname))
  df <- readr::read_csv(fpath, col_types = readr::cols(.default=readr::col_double()))
  
  # last column is outcome
  Y <- df[[ncol(df)]]
  X <- df[, -ncol(df), drop=FALSE]
  
  dataset_label <- sub(".csv$", "", fname)
  
  cat("\n\n=== Running on dataset:", dataset_label, "===\n\n")
  res <- run_comparison_regression(
    dataset_name   = dataset_label,
    X              = X,   # pass data.frame
    Y              = Y,
    drboost_params = DRBOOST_PARAMS
  )
  
  all_results[[dataset_label]] <- res
}

# Combine summary
df_results <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    dataset     = x$dataset,
    DRB_MSE     = x$DRB_MSE,
    DRB_rules   = x$DRB_rules,
    CART_MSE    = x$CART_MSE,
    CART_leaves = x$CART_leaves,
    RF_MSE      = x$RF_MSE,
    GBM_MSE     = x$GBM_MSE,
    XGB_MSE     = x$XGB_MSE,
    stringsAsFactors=FALSE
  )
}))

cat("\nFinal Results:\n")
print(df_results)
cat("\nAll done.\n")
