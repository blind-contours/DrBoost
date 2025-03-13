#' @title Summarize a Disjoint Rule Boosting model
#' @description S3 method for class "disjoint_rule_boost" that prints out
#'   each bounding box, coverage, local model type, and (optionally) 
#'   the mean outcome in that region.
#'
#' @param object   an object of class "disjoint_rule_boost"
#' @param X_train  optional matrix/data.frame of training features 
#'                 (same shape used for disjoint_rule_boost).
#' @param Y_train  optional numeric or binary outcomes (length=nrow(X_train)).
#' @param digits   how many decimal places to show for coefficients
#' @param ...      not used
#'
#' @return invisibly returns a data frame summarizing each rule:
#'   - region_index,
#'   - coverage (# training rows that fall in that bounding box, if X_train is given),
#'   - leftover (TRUE/FALSE if infinite bounding box),
#'   - model_type,
#'   - coefs (short string),
#'   - mean_in_region (if Y_train is provided),
#'   - rule_text (if you used the 'rule_text' in region_info).
#' 
#' @export
summary.disjoint_rule_boost <- function(object,
                                        X_train=NULL,
                                        Y_train=NULL,
                                        digits=3,
                                        ...)
{
  rule_list <- object$rule_list
  if(length(rule_list) == 0) {
    cat("No rules found in this model.\n")
    return(invisible(NULL))
  }
  
  outlist <- list()
  for(i in seq_along(rule_list)) {
    rr <- rule_list[[i]]
    lb <- rr$region_info$lowerBound
    ub <- rr$region_info$upperBound
    mod_type <- rr$region_info$model_type
    
    # check if rule_text exists
    rule_txt <- if(!is.null(rr$region_info$rule_text)) {
      rr$region_info$rule_text
    } else {
      "(no rule_text)"
    }
    
    coefs  <- rr$beta
    # detect leftover box if all infinite
    leftover_flag <- all(is.infinite(lb)) && all(is.infinite(ub))
    
    coverage_i <- NA_integer_
    meanY_i    <- NA_real_
    if(!is.null(X_train)) {
      in_region <- apply(X_train, 1, function(x) {
        all(x >= lb & x <= ub)
      })
      coverage_i <- sum(in_region)
      if(!is.null(Y_train) && coverage_i > 0) {
        meanY_i <- mean(Y_train[in_region])
      }
    }
    
    # Build short string for coefs
    if(length(coefs) == 1) {
      coefs_str <- format(round(coefs, digits), nsmall=digits)
    } else {
      if(length(coefs) > 6) {
        short_vals <- paste(round(coefs[1:6], digits), collapse=", ")
        coefs_str <- paste0(short_vals, ", ...")
      } else {
        coefs_str <- paste(round(coefs, digits), collapse=", ")
      }
    }
    
    outlist[[i]] <- data.frame(
      region_index   = i,
      coverage       = coverage_i,
      leftover       = leftover_flag,
      model_type     = mod_type,
      coefs          = coefs_str,
      mean_in_region = if(is.na(meanY_i)) NA_real_ else round(meanY_i, digits),
      rule_text      = rule_txt,
      stringsAsFactors=FALSE
    )
  }
  
  summary_df <- do.call(rbind, outlist)
  rownames(summary_df) <- NULL
  
  cat("Summary of Disjoint Rule Boosting Model:\n\n")
  print(summary_df)
  
  invisible(summary_df)
}