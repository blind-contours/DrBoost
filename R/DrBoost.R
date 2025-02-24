#' @title Disjoint Rule Boosting (Gradient-Boosting Style) with Excluding Repeated Failing Regions
#' 
#' @description
#'   Builds a piecewise-constant model as a sum of disjoint rules, 
#'   with gradient-boosting style updates and a small loop each iteration 
#'   to skip repeated failing regions. If the top region from beam-search 
#'   hurts validation loss, we exclude that region (so it won't keep appearing)
#'   and re-run beam-search to get the next-best box, up to `max_region_retries` tries.
#'
#' @param X_train matrix or data.frame of shape (n x p)
#' @param Y_train numeric vector (length n): continuous or {0,1}
#' @param X_val matrix or data.frame of shape (m x p)
#' @param Y_val numeric vector (length m)
#' @param outcome_type "continuous" or "binary"
#' @param thresholds_list numeric matrix for bounding-box approach (if logic=FALSE)
#' @param logic boolean; if TRUE, do logic-based search. If FALSE, bounding-box
#' @param K1,K2,K3,K4 integers for beam-search expansions (bounding-box)
#' @param featureNames character vector of length p (for bounding-box rule strings)
#' @param K_twoWay,K_threeWay integers for logic expansions
#' @param max_iterations integer, maximum boosting steps
#' @param learning_rate numeric, e.g. 0.1
#' @param min_obs_pct numeric, min coverage fraction
#' @param max_obs_frac numeric, max coverage fraction to avoid giant rules
#' @param patience integer, early-stopping if no improvement for `patience` steps
#' @param second_order_logistic bool: if TRUE, uses 2nd-order updates for logistic
#' @param max_region_retries integer, how many candidate boxes we try per iteration 
#'        before concluding none helps validation
#'
#' @return A list with:
#'   \item{init}{the initial intercept (numeric scalar)}
#'   \item{rule_list}{list of discovered rules (each with coverage, beta, region info)}
#'   \item{f_train}{final fitted values on training (numeric vector)}
#'   \item{f_val}{final fitted values on validation}
#'   \item{train_loss_final}{final train loss}
#'   \item{val_loss_final}{final validation loss}
#'   \item{validation_losses}{history of validation losses}
#'   \item{best_val_loss}{lowest validation loss found}
#'   \item{best_iteration}{iteration index with best validation}
#'
#' @export
disjoint_rule_boost <- function(
    X_train, Y_train,
    X_val,   Y_val,
    outcome_type = c("continuous","binary"),
    thresholds_list = NULL,
    logic      = FALSE,
    K1         = 200,
    K2         = 200,
    K3         = 200,
    K4         = 200,
    featureNames = character(0),
    K_twoWay   = 100,
    K_threeWay = 100,
    max_iterations = 20,
    learning_rate  = 0.1,
    min_obs_pct    = 0.05,
    max_obs_frac   = 1.0,
    patience       = 3,
    second_order_logistic = TRUE,
    max_region_retries    = 5  # newly added parameter
)
{
  outcome_type <- match.arg(outcome_type, c("continuous","binary"))
  n  <- nrow(X_train)
  m  <- nrow(X_val)
  p  <- ncol(X_train)
  
  X_train_mat <- as.matrix(X_train)
  X_val_mat   <- as.matrix(X_val)
  
  #------------------------------------------
  # 1) Initialize the model (f_train, f_val)
  #------------------------------------------
  f_train <- numeric(n)
  f_val   <- numeric(m)
  
  if (outcome_type == "continuous") {
    # Intercept = mean(Y_train)
    init <- mean(Y_train)
    f_train[] <- init
    f_val[]   <- init
  } else {
    # Intercept = logit(mean(Y_train)) with slight clamp
    eps <- 1e-6
    py  <- mean(Y_train)
    py  <- max(eps, min(1-eps, py))
    init <- log(py / (1 - py))
    f_train[] <- init
    f_val[]   <- init
  }
  
  # Track which points are assigned => no future coverage
  assignedMask <- rep(FALSE, n)
  
  # We'll store discovered rules and their coefficients
  rule_list <- list()
  
  #------------------------------------------
  # 2) Validation metric
  #------------------------------------------
  val_loss_func <- function(y, f_score) {
    if (outcome_type=="continuous") {
      # MSE
      mean((y - f_score)^2)
    } else {
      # logistic => interpret f_score as log-odds
      eps2 <- 1e-15
      p   <- 1 / (1 + exp(-f_score))
      p   <- pmin(pmax(p, eps2), 1 - eps2)
      -mean(y * log(p) + (1-y)*log(1-p))
    }
  }
  
  best_val_loss   <- val_loss_func(Y_val, f_val)
  best_iteration  <- 0
  patience_count  <- 0
  validation_losses <- numeric(0)
  
  #------------------------------------------
  # 3) Iteration Loop
  #------------------------------------------
  for (iter in seq_len(max_iterations)) {
    
    # 3a) Compute negative gradient residual (or 2nd-order approach).
    residuals <- numeric(n)
    if (outcome_type=="continuous") {
      # residual = y - f
      residuals <- (Y_train - f_train)
    } else {
      # logistic
      p_train <- 1 / (1 + exp(-f_train))
      if (!second_order_logistic) {
        # first-order: negative gradient = (y - p)
        residuals <- (Y_train - p_train)
      } else {
        # second-order: store (y - p), then do gamma = sum(y-p)/sum(p(1-p)) after region is found
        residuals <- (Y_train - p_train)
      }
    }
    
    # We'll keep a small loop that tries multiple top-scoring boxes, skipping any that fail.
    excludedBoxes <- list()
    accepted_rule <- FALSE
    local_try <- 1
    
    while (!accepted_rule && (local_try <= max_region_retries)) {
      
      # 3b) Find best region, skipping any that previously failed
      best_reg <- NULL
      
      if (!logic) {
        if (is.null(thresholds_list)) {
          stop("logic=FALSE requires 'thresholds_list'")
        }
        # We call a modified version that can skip excludedBoxes. 
        # Or if your C++ function has not been updated, you have to do it differently.
        best_reg <- find_best_region_beam_cpp(
          X             = X_train_mat,
          residuals     = residuals,
          assignedMask  = assignedMask,
          thresholds_list = thresholds_list,
          min_obs_pct   = min_obs_pct,
          max_obs_frac  = max_obs_frac,
          K1            = K1,
          K2            = K2,
          K3            = K3,
          K4            = K4,
          featureNames  = featureNames,
          excludedBoxes = excludedBoxes # <--- pass excluded region(s)
        )
      } else {
        best_reg <- find_best_binary_rule_3way_topK(
          X          = X_train_mat,
          residuals  = residuals,
          unassigned = which(!assignedMask) - 1L,
          min_obs_pct   = min_obs_pct,
          max_obs_frac  = max_obs_frac,
          K_twoWay      = K_twoWay,
          K_threeWay    = K_threeWay,
          excludedRules = excludedBoxes # <--- similar param for logic-based
        )
      }
      
      if (!isTRUE(best_reg$found_region)) {
        message("No valid region found on attempt ", local_try, " at iteration ", iter, ".")
        break
      }
      
      inside_idx <- best_reg$points_inside + 1L
      coverage   <- length(inside_idx)
      if (coverage == 0) {
        message("Chosen region has coverage=0. Skipping.")
        # exclude it from future tries
        excludedBoxes <- c(excludedBoxes, list(best_reg))
        local_try <- local_try + 1
        next
      }
      
      # 3c) Solve for the coefficient in that region
      #     If continuous => mean of residual in region
      #     If logistic => sum(y-p)/sum(p(1-p)) if second_order_logistic
      new_beta <- 0
      if (outcome_type=="continuous") {
        r_in <- residuals[ inside_idx ]
        new_beta <- mean(r_in)
      } else {
        p_train <- 1 / (1 + exp(-f_train))
        if (!second_order_logistic) {
          r_in <- residuals[ inside_idx ]
          new_beta <- mean(r_in)
        } else {
          eps2 <- 1e-12
          sum_num <- sum(Y_train[inside_idx] - p_train[inside_idx])
          sum_den <- sum(p_train[inside_idx]*(1 - p_train[inside_idx])) + eps2
          new_beta <- sum_num / sum_den
        }
      }
      
      # 3d) Multiply by learning_rate
      new_beta <- learning_rate * new_beta
      
      # 3e) Temporarily update predictions
      old_f_train <- f_train
      old_f_val   <- f_val
      
      # update train
      f_train[ inside_idx ] <- f_train[ inside_idx ] + new_beta
      
      # update val
      if (!logic) {
        lb <- best_reg$lowerBound
        ub <- best_reg$upperBound
        
        for (iv in seq_len(m)) {
          inReg <- TRUE
          for (d in seq_len(p)) {
            vald <- X_val_mat[iv,d]
            if (vald < lb[d] || vald > ub[d]) {
              inReg <- FALSE
              break
            }
          }
          if (inReg) {
            f_val[iv] <- f_val[iv] + new_beta
          }
        }
      } else {
        inside_val_idx <- computeCoverageValLogic(X_val_mat, best_reg$rule_string)
        f_val[ inside_val_idx ] <- f_val[ inside_val_idx ] + new_beta
      }
      
      # 3f) Check validation loss
      val_loss_new <- val_loss_func(Y_val, f_val)
      validation_losses <- c(validation_losses, val_loss_new)
      
      if (val_loss_new < best_val_loss) {
        # accept the rule
        best_val_loss <- val_loss_new
        best_iteration <- iter
        patience_count <- 0
        
        assignedMask[ inside_idx ] <- TRUE
        
        # store the rule
        rule_list[[ length(rule_list)+1 ]] <- list(
          iteration   = iter,
          coverage    = coverage,
          beta        = new_beta,
          region_info = best_reg,
          val_loss    = val_loss_new
        )
        
        cat(sprintf(
          "Iter %2d => coverage=%d, val_loss=%.4f, best_val=%.4f\n  rule: %s\n",
          iter, coverage, val_loss_new, best_val_loss, best_reg$rule_string
        ))
        
        accepted_rule <- TRUE  # stop trying further boxes
      } else {
        # revert
        cat(sprintf(
          "Iter %2d => coverage=%d, val_loss=%.4f (worse), exclude region, revert\n",
          iter, coverage, val_loss_new
        ))
        f_train <- old_f_train
        f_val   <- old_f_val
        # We'll add this failing region to `excludedBoxes`, so next time 
        # the beam-search won't pick it again:
        excludedBoxes <- c(excludedBoxes, list(best_reg))
        
        patience_count <- patience_count + 1
        if (patience_count >= patience) {
          message("Early stopping (no improvement for ", patience, " steps).")
          break
        }
      }
      
      local_try <- local_try + 1
    } # end while over max_region_retries
    
    if (!accepted_rule) {
      # means we didn't accept anything in this iteration 
      message("No acceptable region found at iteration ", iter, ". Stopping.")
      break
    }
  } # end main for(iter)
  
  # final stats
  final_val_loss <- val_loss_func(Y_val, f_val)
  train_loss_final <- if (outcome_type=="continuous") {
    mean((Y_train - f_train)^2)
  } else {
    val_loss_func(Y_train, f_train)
  }
  
  list(
    init              = init,
    rule_list         = rule_list,
    f_train           = f_train,
    f_val             = f_val,
    train_loss_final  = train_loss_final,
    val_loss_final    = final_val_loss,
    validation_losses = validation_losses,
    best_val_loss     = best_val_loss,
    best_iteration    = best_iteration
  )
}
