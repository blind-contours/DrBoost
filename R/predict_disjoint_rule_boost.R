#' @title Predict function for new Disjoint Rule Boosting objects
#'
#' @description
#'   Given the output of disjoint_rule_boost() (which returns
#'   \code{init} + \code{rule_list}), we apply those rules to
#'   new data X_new to compute predictions.
#'
#'   - For outcome_type="continuous", the result is numeric predictions.
#'   - For outcome_type="binary", the result is log-odds by default;
#'     use \code{type="prob"} to get probabilities, or
#'     \code{type="response"} for 0/1 classification.
#'
#' @param object the list returned by disjoint_rule_boost()
#' @param X_new matrix or data.frame of new feature values
#' @param outcome_type either "continuous" or "binary"
#'   (we might store it in \code{object$outcome_type} if you like)
#' @param logic whether the rules are bounding-box or logic
#'   (we could auto-detect if \code{rule_info$region_info$rule_string} is present)
#' @param type "link" (log-odds for binary, numeric for continuous),
#'   "prob" => probability for binary,
#'   "response" => 0/1 for binary.
#'
#' @return Numeric vector of predictions (length = nrow(X_new)).
#'
predict_disjoint_rule_boost <- function(object,
                                        X_new,
                                        outcome_type=c("continuous","binary"),
                                        logic=FALSE,
                                        type=c("link","prob","response")) {
  outcome_type <- match.arg(outcome_type, c("continuous","binary"))
  type         <- match.arg(type, c("link","prob","response"))
  
  # Convert X_new to matrix for bounding-box coverage checks
  X_mat <- as.matrix(X_new)
  ntest <- nrow(X_mat)
  p     <- ncol(X_mat)
  
  # Start prediction at intercept
  out <- rep(object$init, ntest)
  
  # Loop over each discovered rule
  for (rule_obj in object$rule_list) {
    beta_k   <- rule_obj$beta
    reg_info <- rule_obj$region_info
    
    if (!logic) {
      # bounding-box approach
      lb <- reg_info$lowerBound
      ub <- reg_info$upperBound
      
      # Identify which rows in X_new are in [lb[d], ub[d]] for all d
      in_region <- rep(TRUE, ntest)
      for (d in seq_len(p)) {
        xd <- X_mat[, d]
        in_region <- in_region & (xd >= lb[d]) & (xd <= ub[d])
      }
      out[in_region] <- out[in_region] + beta_k
    } else {
      # logic approach => parse rule_string or do "points_inside" for new data
      rule_str <- reg_info$rule_string
      # You need a function that checks coverage of rule_str in X_new
      covered_idx <- computeCoverageValLogic(X_mat, rule_str)
      out[ covered_idx ] <- out[ covered_idx ] + beta_k
    }
  }
  
  if (outcome_type == "continuous") {
    # For regression, out is the final numeric predictions
    if (type != "link") {
      warning("Ignoring 'type', returning raw numeric predictions for continuous outcome.")
    }
    return(out)
  } else {
    # For binary classification, out is the log-odds
    if (type=="link") {
      return(out)  # log-odds
    } else {
      # Convert to probability
      eps <- 1e-15
      p   <- 1 / (1 + exp(-out))
      p   <- pmin(pmax(p, eps), 1 - eps)
      if (type=="prob") {
        return(p)
      } else {
        # type="response" => 0/1
        return(ifelse(p>0.5, 1, 0))
      }
    }
  }
}
