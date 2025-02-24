#!/usr/bin/env Rscript

# ====================================================
# Example R Script: Disjoint Rule Boosting on
# Fingerprinted Chemical Data (Binary), plus optional
# Klekota–Roth JSON dictionary for substructure expansions.
#
# 1) We read a CSV with "KRFP###" columns for fingerprints
#    plus "BIN_XXX" outcome columns.
# 2) We filter fingerprint columns by min_prevalence,
#    rename them X1..Xp for the internal logic-based rules,
#    but keep a rename_map to revert them to KRFP### in output.
# 3) We run disjoint_rule_boost(..., logic=TRUE) for each
#    outcome, storing the final models.
# 4) We produce an extended summary CSV with:
#    - # rules, final val error,
#    - top N rules (by |coef|) with coverage & original KRFP names,
# 5) Optionally, we load a JSON dict that maps each "KRFP###"
#    to a SMARTS pattern, embedding it in the final rule text.
# ====================================================

library(tidyverse)
library(here)
library(devtools)  # if you need load_all() for your local package
load_all()         # loads your local package that defines `disjoint_rule_boost()`
library(stringr)
library(jsonlite)

# ----------------------------------------------------
# 0) Utilities: fix rule strings + optionally embed SMARTS
# ----------------------------------------------------

# (A) Convert "X2" => "KRFP3455" using rename_map, e.g. rename_map["X2"]="KRFP3455"
fixup_rule_string <- function(rule_str, rename_map) {
  if (!nzchar(rule_str)) return(rule_str)
  # Find all occurrences of "X" followed by digits
  matches <- str_extract_all(rule_str, "X\\d+")[[1]]
  if (length(matches) == 0) return(rule_str)
  
  out_str <- rule_str
  for (m in unique(matches)) {
    if (m %in% names(rename_map)) {
      out_str <- gsub(m, rename_map[m], out_str, fixed = TRUE)
    }
  }
  out_str
}

# (B) Convert "KRFP3455" => "KRFP3455=[SMARTS pattern]", if we choose to embed SMARTS
fixup_substructure_string <- function(rule_str, krfp_dict) {
  if (!nzchar(rule_str)) return(rule_str)
  # Find occurrences of "KRFP" followed by digits, e.g. KRFP3752
  matches <- stringr::str_extract_all(rule_str, "KRFP\\d+")[[1]]
  if (length(matches) == 0) return(rule_str)
  
  out_str <- rule_str
  for (m in unique(matches)) {
    if (m %in% names(krfp_dict)) {
      # e.g. substructure could be: "[!#1][CH]([!#1])[!#1]"
      substructure <- krfp_dict[[m]]
      # Instead of "KRFP####=SMILES", we want to replace "KRFP####" with just the substructure
      out_str <- gsub(m, substructure, out_str, fixed = TRUE)
    }
  }
  out_str
}


# ----------------------------------------------------
# 1) Data Loading + Preprocessing
# ----------------------------------------------------

process_data <- function(file_path, min_prevalence = 0.01) {
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Identify fingerprint columns => "KRFP###"
  fp_cols <- grep("^KRFP", colnames(data), value = TRUE)
  X <- as.matrix(data[, fp_cols, drop = FALSE])
  
  # Filter columns by prevalence
  col_prevalence <- colMeans(X == 1, na.rm = TRUE)
  keep_idx <- which(col_prevalence >= min_prevalence)
  if (length(keep_idx) == 0) {
    warning("No fingerprint columns met min_prevalence=", min_prevalence,
            ". Keeping all columns for now.")
  } else if (length(keep_idx) < ncol(X)) {
    message("Dropping ", (ncol(X) - length(keep_idx)),
            " columns with prevalence < ", min_prevalence)
    X <- X[, keep_idx, drop = FALSE]
    fp_cols <- fp_cols[keep_idx]
  }
  
  # Rename columns => X1..Xp
  p_final <- ncol(X)
  new_names <- paste0("X", seq_len(p_final))
  colnames(X) <- new_names
  
  # Build rename_map => "X1" => "KRFP14", etc.
  rename_map <- setNames(fp_cols, new_names)
  
  # Identify outcome columns
  outcome_patterns <- c(
    "BIN_PR","BIN_PXR","BIN_RXR","BIN_GR",
    "BIN_AR","BIN_ERA","BIN_ERB","BIN_FXR",
    "BIN_PPARD","BIN_PPARG","BIN_PPARA"
  )
  
  # "act." => 1, "inact." => 0, else NA
  convert_to_binary <- function(vec) {
    ifelse(vec == "act.", 1, ifelse(vec == "inact.", 0, NA))
  }
  
  outcomes <- list()
  for (pat in outcome_patterns) {
    if (pat %in% colnames(data)) {
      outcomes[[pat]] <- convert_to_binary(data[[pat]])
    } else {
      outcomes[[pat]] <- rep(NA, nrow(data))
      warning("Outcome column not found: ", pat)
    }
  }
  
  list(
    X = X,
    outcomes = outcomes,
    outcome_names = outcome_patterns,
    rename_map = rename_map
  )
}

# ----------------------------------------------------
# 2) Single-Outcome DRB Fit
# ----------------------------------------------------

run_drboost_for_outcome <- function(X, y, receptor_name,
                                    max_iter = 20,
                                    K_twoWay = 50,
                                    K_threeWay = 50,
                                    train_fraction = 0.7,
                                    seed = 123) {
  
  valid_idx <- !is.na(y)
  X_valid <- X[valid_idx, , drop = FALSE]
  y_valid <- y[valid_idx]
  
  if (length(y_valid) < 2) {
    warning("Not enough valid data for receptor=", receptor_name)
    return(NULL)
  }
  set.seed(seed)
  n <- nrow(X_valid)
  train_indices <- sample(seq_len(n), size = floor(train_fraction * n))
  X_train <- X_valid[train_indices, , drop = FALSE]
  y_train <- y_valid[train_indices]
  X_val   <- X_valid[-train_indices, , drop = FALSE]
  y_val   <- y_valid[-train_indices]
  
  # Call disjoint_rule_boost in logic=TRUE mode
  res <- disjoint_rule_boost(
    X_train       = X_train,
    Y_train       = y_train,
    X_val         = X_val,
    Y_val         = y_val,
    outcome_type  = "binary",
    logic         = TRUE,
    max_iterations= max_iter,
    K_twoWay      = K_twoWay,
    K_threeWay    = K_threeWay
  )
  
  list(
    model    = res,
    receptor = receptor_name
  )
}

# ----------------------------------------------------
# 3) Master function: run DRB for each outcome
# ----------------------------------------------------

run_all_receptors <- function(file_path,
                              output_dir     = "drb_results",
                              max_iter       = 20,
                              K_twoWay       = 50,
                              K_threeWay     = 50,
                              min_prevalence = 0.01) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  cat("Loading data from:", file_path, "\n")
  
  dlist     <- process_data(file_path, min_prevalence = min_prevalence)
  X         <- dlist$X
  outcomes  <- dlist$outcomes
  out_names <- dlist$outcome_names
  rename_map<- dlist$rename_map
  
  results <- list()
  for (receptor in out_names) {
    cat("\n[INFO] Fitting DRB for outcome:", receptor, "\n")
    y_vec <- outcomes[[receptor]]
    
    fit_result <- run_drboost_for_outcome(
      X = X, y = y_vec, receptor_name = receptor,
      max_iter = max_iter, K_twoWay = K_twoWay, K_threeWay = K_threeWay
    )
    if (is.null(fit_result)) {
      next
    }
    # Store final model + rename_map
    results[[receptor]] <- list(
      model      = fit_result$model,
      receptor   = receptor,
      rename_map = rename_map
    )
    
    # Save each outcome model
    rds_path <- file.path(output_dir, paste0(receptor, "_model.rds"))
    saveRDS(results[[receptor]], file = rds_path)
    
    # Optionally write discovered rules
    if (!is.null(fit_result$model$rule_list)) {
      discovered_rules <- fit_result$model$rule_list
      rule_summary <- data.frame(
        iteration   = seq_along(discovered_rules),
        rule_string = sapply(discovered_rules, function(rr) rr$rule_string),
        coverage    = sapply(discovered_rules, function(rr)
          if (!is.null(rr$coverage)) sum(rr$coverage) else NA),
        coef        = sapply(discovered_rules, function(rr)
          if (!is.null(rr$coef)) rr$coef else NA),
        stringsAsFactors = FALSE
      )
      csv_path <- file.path(output_dir, paste0(receptor, "_rules.csv"))
      write.csv(rule_summary, file = csv_path, row.names = FALSE)
    }
  }
  
  # Save global results
  saveRDS(results, file = file.path(output_dir, "all_models.rds"))
  cat("\nAll receptor models finished. Results saved in:", output_dir, "\n")
  results
}

# ----------------------------------------------------
# 4) Summarize with top N rules, optionally mapping
#    KRFP### => SMARTS from a JSON dictionary
# ----------------------------------------------------

summarize_drb_models_extended <- function(results,
                                          output_dir   = "drb_results",
                                          n_top        = 3,
                                          krfp_json    = NULL,
                                          use_smarts   = FALSE) {
  # If provided, load the KRFP dictionary (JSON) => e.g. "KRFP14" => "[!#1][CH]..."
  krfp_dict <- NULL
  if (!is.null(krfp_json) && file.exists(krfp_json)) {
    cat("Loading KRFP dictionary from:", krfp_json, "\n")
    krfp_dict <- fromJSON(krfp_json)  # e.g. list(KRFP1="[!#1][CH]..." etc.)
  }
  
  summary_list <- list()
  
  for (receptor in names(results)) {
    res_obj <- results[[receptor]]
    if (is.null(res_obj$model)) next
    
    model_fit  <- res_obj$model
    rename_map <- res_obj$rename_map
    rule_list  <- model_fit$rule_list
    
    final_coefs <- coef(model_fit$final_model)
    if ("Intercept_Col" %in% names(final_coefs)) {
      final_coefs <- final_coefs[ setdiff(names(final_coefs), "Intercept_Col") ]
    }
    
    n_rules      <- if (!is.null(rule_list)) length(rule_list) else 0
    final_val_err<- if (!is.null(model_fit$val_mse)) model_fit$val_mse else NA_real_
    
    # Sort by abs value, descending
    sorted_cols <- names(final_coefs)[order(abs(final_coefs), decreasing = TRUE)]
    n_show <- min(length(sorted_cols), n_top)
    
    out_row <- list(
      receptor       = receptor,
      num_rules      = n_rules,
      final_val_loss = final_val_err
    )
    
    # Match coefficient => rule coverage => fix up "X##" => "KRFP###"
    # then optionally => SMARTS
    get_rule_info <- function(coln) {
      if (is.null(rule_list)) return(list(rule=coln, cover=NA))
      coln_clean <- gsub("`", "", coln)
      wh <- which(sapply(rule_list, function(rr) rr$rule_string) == coln_clean)
      if (length(wh)>0) {
        idx       <- wh[1]
        raw_rule  <- rule_list[[idx]]$rule_string
        coverage  <- rule_list[[idx]]$coverage
        # Step1: "X##" => "KRFP###"
        rule_fixed <- fixup_rule_string(raw_rule, rename_map)
        # Step2 (optional): "KRFP###" => "KRFP###=[SMARTS]"
        if (use_smarts && !is.null(krfp_dict)) {
          rule_fixed <- fixup_substructure_string(rule_fixed, krfp_dict)
        }
        return(list(rule = rule_fixed, cover = coverage))
      }
      list(rule=coln, cover=NA)
    }
    
    for (i in seq_len(n_show)) {
      coln    <- sorted_cols[i]
      info    <- get_rule_info(coln)
      the_coef<- final_coefs[coln]
      
      out_row[[paste0("topRule", i)]] <- info$rule
      out_row[[paste0("topCoef", i)]] <- the_coef
      out_row[[paste0("topCov",  i)]] <- info$cover
    }
    summary_list[[receptor]] <- out_row
  }
  
  if (length(summary_list)==0) {
    summary_df <- data.frame()
  } else {
    summary_df <- do.call(rbind.data.frame, summary_list)
  }
  out_csv <- file.path(output_dir, "summary_all_receptors.csv")
  write.csv(summary_df, out_csv, row.names = FALSE)
  summary_df
}

# ----------------------------------------------------
# 5) Main usage example
# ----------------------------------------------------

if (sys.nframe() == 0) {
  csv_path <- here("data/Nura_database_merged_fingerprints.csv")
  out_dir  <- "drb_results"
  
  # Hyperparameters
  max_iter      <- 20
  K2            <- 50
  K3            <- 50
  min_prevalent <- 0.1
  
  # 1) Fit models for all receptors
  res_all <- run_all_receptors(
    file_path      = csv_path,
    output_dir     = out_dir,
    max_iter       = max_iter,
    K_twoWay       = K2,
    K_threeWay     = K3,
    min_prevalence = min_prevalent
  )
  
  # 2) Summarize with top 2 rules => 
  #    'X##' => 'KRFP###', and optionally 'KRFP###' => SMARTS if use_smarts=TRUE
  #    Provide a JSON file if you want the full substructure expansions
  my_krfp_json <- here("data", "krfp.json")  # adjust path as needed
  
  final_summary <- summarize_drb_models_extended(
    results     = res_all,
    output_dir  = out_dir,
    n_top       = 1,
    krfp_json   = my_krfp_json,  # can be NULL if you don't have one
    use_smarts  = TRUE          # set to FALSE if you only want short KRFP labels
  )
  print(final_summary)
}
