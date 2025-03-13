# -------------------------------------------------------------------
# Helper to build text rule from lb,ub
# -------------------------------------------------------------------
build_box_string <- function(lb, ub, feat_names, digits=3) {
  p <- length(lb)
  parts <- character(0)
  if(is.null(feat_names) || length(feat_names)!=p) {
    feat_names <- paste0("X", seq_len(p))
  }
  for(j in seq_len(p)) {
    conds_j <- character(0)
    if(!is.infinite(lb[j])) {
      conds_j <- c(conds_j, paste0(
        feat_names[j]," >= ", round(lb[j], digits)
      ))
    }
    if(!is.infinite(ub[j])) {
      conds_j <- c(conds_j, paste0(
        feat_names[j]," <= ", round(ub[j], digits)
      ))
    }
    if(length(conds_j)>0) {
      parts <- c(parts, paste(conds_j, collapse=" & "))
    }
  }
  if(length(parts)==0) {
    return("All features: (-Inf, +Inf)")
  } else {
    return(paste(parts, collapse=" & "))
  }
}

# -------------------------------------------------------------------
# Build local design matrix (constant, linear, poly2)
# -------------------------------------------------------------------
build_local_design <- function(Xsub, local_model_type = "constant") {
  nsub <- nrow(Xsub)
  p    <- ncol(Xsub)
  if (nsub < 1 || p < 1) {
    return(matrix(1, nrow = nsub, ncol = 1))
  }
  if (local_model_type == "constant") {
    return(matrix(1, nrow = nsub, ncol = 1))
  } else if (local_model_type == "linear") {
    return(cbind(1, Xsub))
  } else if (local_model_type == "poly2") {
    out <- cbind(1, Xsub)
    for (j in seq_len(p)) {
      out <- cbind(out, Xsub[, j]^2)
    }
    return(out)
  }
  # fallback
  return(matrix(1, nrow = nsub, ncol = 1))
}

# -------------------------------------------------------------------
# Fit local submodel + measure SSE improvement
# Returns: list(coeffs=..., model_type=..., sse_improve=...)
# -------------------------------------------------------------------
fit_local_submodel_sse <- function(X_in, r_in,
                                   candidate_types=c("constant","linear","poly2"),
                                   base_sse=NULL) {
  if(is.null(base_sse)) {
    base_sse <- sum(r_in^2)
  }
  best_fit <- list(coeffs=numeric(0), model_type=NA_character_, sse_improve=-Inf)
  if(nrow(X_in) < 1) return(best_fit)
  
  for(tt in candidate_types) {
    des <- build_local_design(X_in, tt)
    olsfit <- tryCatch(stats::lm.fit(des, r_in), error=function(e)NULL)
    if(is.null(olsfit) || anyNA(olsfit$coefficients)) next
    
    pred_loc <- as.vector(des %*% olsfit$coefficients)
    new_sse  <- sum((r_in - pred_loc)^2)
    improvement <- base_sse - new_sse
    if(improvement > best_fit$sse_improve) {
      best_fit$sse_improve = improvement
      best_fit$model_type  = tt
      best_fit$coeffs      = olsfit$coefficients
    }
  }
  return(best_fit)
}

# -------------------------------------------------------------------
# Restrict thresholds to be within [parent_lb, parent_ub]
# Then call find_topk_candidates_beam_cpp
# -------------------------------------------------------------------
restrict_thresholds_and_search <- function(
    X_sub, r_sub, coverageCount_sub,
    parent_lb, parent_ub,
    max_overlap,
    thresholds_global,
    min_obs_pct,
    max_obs_frac,
    K1,K2,K3,K4,
    topK
) {
  p <- ncol(X_sub)
  # Make a copy
  tmat <- thresholds_global
  for(j in seq_len(p)) {
    oldcol <- tmat[, j]
    valid  <- (oldcol >= parent_lb[j]) & (oldcol <= parent_ub[j])
    # mark invalid thresholds as NA
    tmat[!valid, j] <- NA_real_
  }
  
  out <- find_topk_candidates_beam_cpp(
    X_sub, r_sub, coverageCount_sub, max_overlap,
    thresholds_list = tmat,
    min_obs_pct = min_obs_pct,
    max_obs_frac= max_obs_frac,
    K1=K1, K2=K2, K3=K3, K4=K4,
    topK=topK
  )
  return(out)
}

# -------------------------------------------------------------------
# pick_best_box_and_model: from a topK list, do SSE improvement
# -------------------------------------------------------------------
pick_best_box_and_model <- function(
    candList, X, r, coverageCount,
    coverage_idx_map = NULL,
    outcome_type=c("continuous","binary"),
    learning_rate=0.1
) {
  outcome_type <- match.arg(outcome_type, c("continuous","binary"))
  if(!candList[["found_any"]]) return(NULL)
  lbMat <- candList[["lowerBounds"]]
  ubMat <- candList[["upperBounds"]]
  insideList <- candList[["inside"]]
  nCands <- nrow(lbMat)
  if(nCands < 1) return(NULL)
  
  best_gain <- -Inf
  best_item <- NULL
  
  for(iCand in seq_len(nCands)) {
    lbVec  <- lbMat[iCand, ]
    ubVec  <- ubMat[iCand, ]
    insideIdx <- insideList[[iCand]]  # 0-based
    if(!is.null(coverage_idx_map)) {
      # map to global
      insideIdx <- coverage_idx_map[ insideIdx + 1 ]
    } else {
      insideIdx <- insideIdx + 1L
    }
    
    if(length(insideIdx)<1) next
    X_in <- X[insideIdx, , drop=FALSE]
    r_in <- r[insideIdx]
    base_sse <- sum(r_in^2)
    
    if(outcome_type=="binary") {
      # piecewise-constant only
      rawBeta <- mean(r_in)
      improvement <- base_sse - sum((r_in - rawBeta)^2)
      if(improvement> best_gain) {
        best_gain <- improvement
        best_item <- list(
          lb=lbVec, ub=ubVec,
          coverage=length(insideIdx),
          model_type="binary-constant",
          coeffs= rawBeta
        )
      }
    } else {
      local_fit <- fit_local_submodel_sse(
        X_in, r_in,
        candidate_types=c("constant","linear","poly2"),
        base_sse = base_sse
      )
      if(local_fit$sse_improve > best_gain) {
        best_gain <- local_fit$sse_improve
        best_item <- list(
          lb=lbVec, ub=ubVec,
          coverage=length(insideIdx),
          model_type= local_fit$model_type,
          coeffs= local_fit$coeffs
        )
      }
    }
  }
  
  if(is.null(best_item)) return(NULL)
  # apply learning rate
  best_item$coeffs <- best_item$coeffs * learning_rate
  best_item
}

# -------------------------------------------------------------------
# A helper to measure the validation loss
# (for continuous or binary)
# -------------------------------------------------------------------
make_loss_function <- function(Y_val, outcome_type) {
  if(outcome_type=="continuous") {
    function(f_val) mean((Y_val - f_val)^2)
  } else {
    eps <- 1e-6
    function(f_val) {
      p_ <- 1/(1+exp(-f_val))
      p_ <- pmin(pmax(p_, eps), 1-eps)
      -mean(Y_val*log(p_) + (1-Y_val)*log(1-p_))
    }
  }
}

# -------------------------------------------------------------------
# Recursive leftover subdivision
# -------------------------------------------------------------------
recursive_subdivide_leftover <- function(
    X, Y,
    f_train,       # current model predictions on train
    coverageCount,
    max_overlap=1,
    thresholds_list,
    outcome_type=c("continuous","binary"),
    learning_rate=0.1,
    min_obs_pct=0.01,
    max_obs_frac=1.0,
    K1=300,K2=300,K3=300,K4=300,
    topK_candidates=50,
    leftover_patience=5,
    max_iter=100,
    val_loss_func=NULL,   # pass a function(Y_val, f_val) -> numeric
    f_val=NULL,           # the current validation predictions
    best_val_loss=NULL
) {
  outcome_type <- match.arg(outcome_type, c("continuous","binary"))
  leftoverBoxes <- list()
  no_improve_count <- 0
  iter <- 1
  n <- nrow(X)
  # We'll define a function to measure training "loss" or SSE if continuous
  measure_train_loss <- NULL
  if(outcome_type=="continuous") {
    measure_train_loss <- function(y, f) sum((y - f)^2)
  } else {
    # logistic deviance
    measure_train_loss <- function(y, f) {
      eps <- 1e-9
      p_ <- 1/(1+exp(-f))
      p_ <- pmin(pmax(p_, eps), 1-eps)
      -sum(y*log(p_) + (1-y)*log(1-p_))
    }
  }
  base_train_loss <- measure_train_loss(Y, f_train)
  if(is.null(best_val_loss)) best_val_loss <- if(!is.null(val_loss_func)) val_loss_func(f_val) else Inf
  
  repeat {
    leftover_idx <- which(coverageCount < max_overlap)
    if(length(leftover_idx)<5) {
      message("Recursive leftover: too few uncovered points => stop.")
      break
    }
    # get residual
    if(outcome_type=="continuous") {
      r_ <- Y - f_train
    } else {
      p_ <- 1/(1+exp(-f_train))
      r_ <- Y - p_
    }
    
    # gather topK
    candList <- find_topk_candidates_beam_cpp(
      X[leftover_idx,,drop=FALSE],
      r_[leftover_idx],
      coverageCount= coverageCount[leftover_idx],
      max_overlap=0,  # or use max_overlap=0 for leftover sub-problem
      thresholds_list= thresholds_list,
      min_obs_pct= min_obs_pct,
      max_obs_frac= max_obs_frac,
      K1=K1, K2=K2, K3=K3, K4=K4,
      topK= topK_candidates
    )
    if(!candList$found_any) {
      message("Recursive leftover: no subregion found => stop.")
      break
    }
    
    chosen_box <- pick_best_box_and_model(
      candList,
      X, r_,
      coverageCount,
      coverage_idx_map = leftover_idx - 1,
      outcome_type= outcome_type,
      learning_rate= learning_rate
    )
    if(is.null(chosen_box)) {
      message("Recursive leftover: no improvement from submodels => break.")
      break
    }
    
    lb <- chosen_box$lb
    ub <- chosen_box$ub
    X_subLeft <- X[leftover_idx,,drop=FALSE]
    inside_sub <- which(apply(X_subLeft,1,function(z) all(z>=lb & z<=ub)))
    if(length(inside_sub)<1) {
      message("Chosen leftover subregion coverage=0 => skip.")
      break
    }
    inside_global <- leftover_idx[ inside_sub ]
    
    old_f_train <- f_train
    old_f_val   <- f_val
    old_val_loss<- best_val_loss
    
    # apply local shift
    if(outcome_type=="binary") {
      shiftVal <- chosen_box$coeffs[1]
      f_train[inside_global] <- f_train[inside_global] + shiftVal
    } else {
      X_in <- X[inside_global,,drop=FALSE]
      des  <- build_local_design(X_in, chosen_box$model_type)
      pred <- as.vector(des %*% chosen_box$coeffs)
      f_train[inside_global] <- f_train[inside_global] + pred
    }
    
    # also update validation if we have f_val
    if(!is.null(f_val)) {
      inside_val <- which(apply(as.matrix(X_val),1,function(z){
        all(z>=lb & z<=ub)
      }))
      if(length(inside_val)>0) {
        if(outcome_type=="binary") {
          shiftVal <- chosen_box$coeffs[1]
          f_val[inside_val] <- f_val[inside_val] + shiftVal
        } else {
          X_val_sub <- X_val[ inside_val, , drop=FALSE ]
          dval <- build_local_design(X_val_sub, chosen_box$model_type)
          pval <- as.vector(dval %*% chosen_box$coeffs)
          f_val[inside_val] <- f_val[inside_val] + pval
        }
      }
    }
    
    # coverage update
    coverageCount[inside_global] <- coverageCount[inside_global] + 1L
    
    # check val improvement
    new_val_loss <- if(!is.null(val_loss_func)) val_loss_func(f_val) else Inf
    
    if(new_val_loss < best_val_loss) {
      best_val_loss <- new_val_loss
      no_improve_count <- 0
      leftoverBoxes[[ length(leftoverBoxes)+1 ]] <- list(
        iteration = iter,
        lowerBound= lb,
        upperBound= ub,
        model_type= chosen_box$model_type,
        coeffs    = as.numeric(chosen_box$coeffs),
        coverage  = length(inside_global)
      )
      message(sprintf("Recursive leftover => coverage=%d => val_loss=%.4f (improved)",
                      length(inside_global), new_val_loss))
    } else {
      message(sprintf("Recursive leftover => coverage=%d => val_loss=%.4f (worse), revert",
                      length(inside_global), new_val_loss))
      # revert
      f_train <- old_f_train
      if(!is.null(f_val)) f_val <- old_f_val
      coverageCount[inside_global] <- coverageCount[inside_global] - 1L
      no_improve_count <- no_improve_count + 1
      if(no_improve_count >= leftover_patience) {
        message("Recursive leftover => early stop (no improvement).")
        break
      }
    }
    
    iter <- iter + 1
    if(iter > max_iter) {
      message("Recursive leftover => reached max iteration => stop.")
      break
    }
  }
  
  list(
    coverageCount   = coverageCount,
    leftoverBoxes   = leftoverBoxes,
    f_train         = f_train,
    f_val           = f_val,
    best_val_loss   = best_val_loss
  )
}

# -------------------------------------------------------------------
# disjoint_rule_boost: pass1, pass2, then recursive leftover
# -------------------------------------------------------------------
disjoint_rule_boost <- function(
    X_train, Y_train,
    X_val,   Y_val,
    outcome_type  = c("continuous","binary"),
    thresholds_list,
    max_pass1_iter=100,
    max_dim_refine=4,
    patience      =10,
    learning_rate =0.1,
    min_obs_pct   =0.01,
    max_obs_frac  =1.0,
    K1=300, K2=300, K3=300, K4=300,
    max_overlap   =1,
    refine_leftover=TRUE,
    topK_candidates=50,
    leftover_patience=5,  # patience for leftover recursion
    ...
) {
  outcome_type <- match.arg(outcome_type, c("continuous","binary"))
  X_train <- as.matrix(X_train)
  X_val   <- as.matrix(X_val)
  n <- nrow(X_train); p <- ncol(X_train)
  m <- nrow(X_val)
  
  feat_names <- colnames(X_train)
  
  # Initialize preds + loss
  if(outcome_type=="continuous") {
    init <- mean(Y_train)
    f_train <- rep(init, n)
    f_val   <- rep(init, m)
  } else {
    eps <- 1e-6
    py  <- mean(Y_train)
    py  <- max(eps, min(1-eps, py))
    init <- log(py/(1-py))
    f_train <- rep(init, n)
    f_val   <- rep(init, m)
  }
  loss_func <- make_loss_function(Y_val, outcome_type)
  best_val_loss <- loss_func(f_val)
  val_loss_history <- numeric(0)
  val_loss_history <- c(val_loss_history, best_val_loss)
  patience_count <- 0
  
  coverageCount <- integer(n)
  region_list   <- list()
  
  # helper for pass1 picks
  pick_best_box_and_model_pass1 <- function(candList, X, r, coverageCount, coverage_idx_map=NULL) {
    pick_best_box_and_model(
      candList, X, r, coverageCount,
      coverage_idx_map = coverage_idx_map,
      outcome_type= outcome_type,
      learning_rate= learning_rate
    )
  }
  
  # ----------------------------------------------------------------
  # PASS 1
  # ----------------------------------------------------------------
  for(iter in seq_len(max_pass1_iter)) {
    if(outcome_type=="continuous") {
      residuals <- Y_train - f_train
    } else {
      p_ <- 1/(1+exp(-f_train))
      residuals <- Y_train - p_
    }
    
    candList <- find_topk_candidates_beam_cpp(
      X_train, residuals,
      coverageCount = coverageCount,
      max_overlap   = max_overlap,
      thresholds_list= thresholds_list,
      min_obs_pct   = min_obs_pct,
      max_obs_frac  = max_obs_frac,
      K1=K1, K2=K2, K3=ifelse(max_dim_refine>=3,K3,0), K4=ifelse(max_dim_refine>=4,K4,0),
      topK=topK_candidates
    )
    if(!candList$found_any) {
      message("Pass1: no valid region found at iteration ", iter)
      break
    }
    
    chosen_box <- pick_best_box_and_model_pass1(
      candList, X_train, residuals, coverageCount
    )
    if(is.null(chosen_box)) {
      message("Pass1: no improvement => break.")
      break
    }
    
    old_f_train <- f_train
    old_f_val   <- f_val
    lb <- chosen_box$lb
    ub <- chosen_box$ub
    
    inside_global <- which(apply(X_train,1,function(z){
      all(z>=lb & z<=ub)
    }))
    if(length(inside_global)<1) {
      message("Pass1 coverage=0 => skip.")
      break
    }
    
    # apply shift
    if(outcome_type=="binary") {
      shiftVal <- chosen_box$coeffs[1]
      f_train[inside_global] <- f_train[inside_global] + shiftVal
    } else {
      X_in <- X_train[inside_global, , drop=FALSE]
      des  <- build_local_design(X_in, chosen_box$model_type)
      pred_loc <- as.vector(des %*% chosen_box$coeffs)
      f_train[inside_global] <- f_train[inside_global] + pred_loc
    }
    
    # update val
    inside_val <- which(apply(X_val,1,function(z){
      all(z>=lb & z<=ub)
    }))
    if(length(inside_val)>0) {
      if(outcome_type=="binary") {
        shiftVal <- chosen_box$coeffs[1]
        f_val[inside_val] <- f_val[inside_val] + shiftVal
      } else {
        X_val_sub <- X_val[inside_val,,drop=FALSE]
        des_val_sub<- build_local_design(X_val_sub, chosen_box$model_type)
        pred_val_sub<- as.vector(des_val_sub %*% chosen_box$coeffs)
        f_val[inside_val] <- f_val[inside_val] + pred_val_sub
      }
    }
    
    new_val_loss <- loss_func(f_val)
    val_loss_history <- c(val_loss_history, new_val_loss)
    
    if(new_val_loss < best_val_loss) {
      best_val_loss <- new_val_loss
      patience_count <- 0
      coverageCount[inside_global] <- coverageCount[inside_global] + 1L
      region_list[[ length(region_list)+1 ]] <- list(
        iteration  = iter,
        coverage   = length(inside_global),
        lowerBound = lb,
        upperBound = ub,
        model_type = chosen_box$model_type,
        coeffs     = as.numeric(chosen_box$coeffs)
      )
      message(sprintf("Pass1 iter %2d => coverage=%d, val=%.4f (improved)",
                      iter, length(inside_global), new_val_loss))
    } else {
      message(sprintf("Pass1 iter %2d => coverage=%d, val=%.4f (worse), revert",
                      iter, length(inside_global), new_val_loss))
      f_train <- old_f_train
      f_val   <- old_f_val
      patience_count <- patience_count + 1
      if(patience_count>=patience) {
        message("Early stop pass1 => no improvement.")
        break
      }
    }
  } # pass1
  
  # ----------------------------------------------------------------
  # PASS 2: refine each box
  # ----------------------------------------------------------------
  refinedBoxes <- list()
  for(iBox in seq_along(region_list)) {
    parentBox <- region_list[[iBox]]
    lbp <- parentBox$lowerBound
    ubp <- parentBox$upperBound
    inside_idx <- which(apply(X_train,1,function(z){
      all(z>=lbp & z<=ubp)
    }))
    if(length(inside_idx)<5) {
      refinedBoxes[[iBox]] <- list(parentBox=parentBox, subBoxes=NULL)
      next
    }
    subBoxes <- list()
    
    for(localIter in seq_len(5)) {
      if(outcome_type=="continuous") {
        r_local <- Y_train - f_train
      } else {
        p_ <- 1/(1+exp(-f_train))
        r_local <- Y_train - p_
      }
      X_sub2   <- X_train[inside_idx,,drop=FALSE]
      r_sub2   <- r_local[inside_idx]
      cCount_sub2 <- coverageCount[inside_idx]
      
      subCand <- restrict_thresholds_and_search(
        X_sub2, r_sub2, cCount_sub2,
        parent_lb=lbp, parent_ub=ubp,
        max_overlap=max_overlap,
        thresholds_global=thresholds_list,
        min_obs_pct=min_obs_pct,
        max_obs_frac=max_obs_frac,
        K1=K1, K2=K2,
        K3=ifelse(max_dim_refine>=3, K3,0),
        K4=ifelse(max_dim_refine>=4, K4,0),
        topK=topK_candidates
      )
      if(!subCand$found_any) {
        message(sprintf("No subregion found => box %d localIter=%d", iBox, localIter))
        break
      }
      chosen_sub <- pick_best_box_and_model(
        subCand, X_train, r_local,
        coverageCount,
        coverage_idx_map = inside_idx - 1,
        outcome_type= outcome_type,
        learning_rate= learning_rate
      )
      if(is.null(chosen_sub)) {
        message(sprintf("No subregion improvement => stop refine (box %d)", iBox))
        break
      }
      
      old_f_train2 <- f_train
      old_f_val2   <- f_val
      sub_lb <- chosen_sub$lb
      sub_ub <- chosen_sub$ub
      
      inside_sub_global <- inside_idx[ which(apply(X_sub2,1,function(z){
        all(z>=sub_lb & z<=sub_ub)
      })) ]
      if(length(inside_sub_global)<1) {
        message("Refine coverage=0 => skip.")
        break
      }
      
      # apply shift
      if(outcome_type=="binary") {
        shiftVal2 <- chosen_sub$coeffs[1]
        f_train[ inside_sub_global ] <- f_train[ inside_sub_global ] + shiftVal2
      } else {
        X_in2 <- X_train[ inside_sub_global, , drop=FALSE]
        des_in2 <- build_local_design(X_in2, chosen_sub$model_type)
        pred2 <- as.vector(des_in2 %*% chosen_sub$coeffs)
        f_train[ inside_sub_global ] <- f_train[ inside_sub_global ] + pred2
      }
      
      # update val
      inside_val_sub <- which(apply(X_val,1,function(z){
        all(z>=sub_lb & z<=sub_ub)
      }))
      if(length(inside_val_sub)>0) {
        if(outcome_type=="binary") {
          shiftVal2 <- chosen_sub$coeffs[1]
          f_val[ inside_val_sub ] <- f_val[ inside_val_sub ] + shiftVal2
        } else {
          X_val2 <- X_val[ inside_val_sub, , drop=FALSE ]
          des_val2 <- build_local_design(X_val2, chosen_sub$model_type)
          pred_val2<- as.vector(des_val2 %*% chosen_sub$coeffs)
          f_val[ inside_val_sub ] <- f_val[ inside_val_sub ] + pred_val2
        }
      }
      
      new_vloss2 <- loss_func(f_val)
      if(new_vloss2 < best_val_loss) {
        best_val_loss <- new_vloss2
        coverageCount[ inside_sub_global ] <-
          coverageCount[ inside_sub_global ] + 1L
        subBoxes[[ length(subBoxes)+1 ]] <- list(
          coverage   = length(inside_sub_global),
          lowerBound = sub_lb,
          upperBound = sub_ub,
          model_type = chosen_sub$model_type,
          coeffs     = as.numeric(chosen_sub$coeffs)
        )
        message(sprintf("Refine => improved val=%.4f => keep (box %d)",
                        new_vloss2, iBox))
      } else {
        message(sprintf("Refine => val=%.4f => revert (box %d)",
                        new_vloss2, iBox))
        f_train <- old_f_train2
        f_val   <- old_f_val2
        break
      }
    }
    refinedBoxes[[iBox]] <- list(parentBox=parentBox, subBoxes=subBoxes)
  }
  
  # combine pass1 + pass2 into a preliminary rule_list
  rule_list <- list()
  for(rb in region_list) {
    if(is.null(rb$lowerBound)) next
    rule_list[[ length(rule_list)+1 ]] <- list(
      region_info = list(
        lowerBound= rb$lowerBound,
        upperBound= rb$upperBound,
        model_type= rb$model_type,
        rule_text = build_box_string(rb$lowerBound, rb$upperBound, feat_names)
      ),
      beta = rb$coeffs
    )
  }
  # refined sub-boxes
  for(xrf in refinedBoxes) {
    if(length(xrf$subBoxes)<1) next
    for(sb in xrf$subBoxes) {
      rule_list[[ length(rule_list)+1 ]] <- list(
        region_info= list(
          lowerBound= sb$lowerBound,
          upperBound= sb$upperBound,
          model_type= sb$model_type,
          rule_text = build_box_string(sb$lowerBound, sb$upperBound, feat_names)
        ),
        beta= sb$coeffs
      )
    }
  }
  
  # ----------------------------------------------------------------
  # RECURSIVE LEFTOVER SUBDIVISION
  # ----------------------------------------------------------------
  leftoverBoxes <- NULL
  if(refine_leftover) {
    resLeft <- recursive_subdivide_leftover(
      X = X_train, Y = Y_train,
      f_train= f_train,
      coverageCount= coverageCount,
      max_overlap = max_overlap,
      thresholds_list= thresholds_list,
      outcome_type= outcome_type,
      learning_rate= learning_rate,
      min_obs_pct= min_obs_pct,
      max_obs_frac= max_obs_frac,
      K1=K1, K2=K2, K3=K3, K4=K4,
      topK_candidates= topK_candidates,
      leftover_patience= leftover_patience,
      max_iter= 100,
      val_loss_func= loss_func,
      f_val= f_val,
      best_val_loss= best_val_loss
    )
    coverageCount <- resLeft$coverageCount
    leftoverBoxes <- resLeft$leftoverBoxes
    f_train       <- resLeft$f_train
    f_val         <- resLeft$f_val
    best_val_loss <- resLeft$best_val_loss
    
    # add leftoverBoxes to rule_list
    if(length(leftoverBoxes)>0) {
      for(lbobj in leftoverBoxes) {
        rule_list[[ length(rule_list)+1 ]] <- list(
          region_info= list(
            lowerBound= lbobj$lowerBound,
            upperBound= lbobj$upperBound,
            model_type= lbobj$model_type,
            rule_text = build_box_string(lbobj$lowerBound, lbobj$upperBound, feat_names)
          ),
          beta= as.numeric(lbobj$coeffs)
        )
      }
    }
  }
  
  # final
  val_loss_history <- c(val_loss_history, best_val_loss)
  out_obj <- list(
    init = init,
    coverageCount   = coverageCount,
    region_list     = region_list,
    leftoverBoxes   = leftoverBoxes,
    refinedBoxes    = refinedBoxes,
    rule_list       = rule_list,
    f_train_final   = f_train,
    f_val_final     = f_val,
    val_loss_history= val_loss_history,
    best_val_loss   = best_val_loss
  )
  class(out_obj) <- "disjoint_rule_boost"
  return(out_obj)
}

# -------------------------------------------------------------------
# Predict function for disjoint_rule_boost
# -------------------------------------------------------------------
predict_disjoint_rule_boost <- function(
    object, X_new,
    outcome_type=c("continuous","binary"),
    type=c("link","prob","response")
) {
  outcome_type <- match.arg(outcome_type, c("continuous","binary"))
  type         <- match.arg(type, c("link","prob","response"))
  X_mat <- as.matrix(X_new)
  ntest <- nrow(X_mat)
  p     <- ncol(X_mat)
  pred  <- rep(object$init, ntest)
  
  for(rule_obj in object$rule_list) {
    beta_vec  <- rule_obj$beta
    reg_info  <- rule_obj$region_info
    lb        <- reg_info$lowerBound
    ub        <- reg_info$upperBound
    mod_type  <- reg_info$model_type
    
    in_region <- rep(TRUE, ntest)
    for(d in seq_len(p)) {
      xd <- X_mat[, d]
      in_region <- in_region & (xd >= lb[d]) & (xd <= ub[d])
    }
    if(!any(in_region)) next
    
    if(length(beta_vec)==1) {
      pred[in_region] <- pred[in_region] + beta_vec
    } else {
      X_sub <- X_mat[in_region, , drop=FALSE]
      des   <- build_local_design(X_sub, mod_type)
      local_pred <- as.vector(des %*% beta_vec)
      pred[in_region] <- pred[in_region] + local_pred
    }
  }
  
  if(outcome_type=="continuous") {
    if(type!="link") {
      warning("For continuous outcome, returning numeric predictions; 'type' ignored.")
    }
    return(pred)
  } else {
    # outcome_type=="binary" => log-odds
    if(type=="link") {
      return(pred)
    } else {
      eps  <- 1e-15
      pval <- 1/(1+exp(-pred))
      pval <- pmin(pmax(pval, eps), 1-eps)
      if(type=="prob") {
        return(pval)
      } else {
        return(ifelse(pval>0.5, 1, 0))
      }
    }
  }
}
