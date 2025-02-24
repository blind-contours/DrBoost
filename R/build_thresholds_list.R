#' @title Build Thresholds List for Disjoint Rule Boosting
#'
#' @description
#'   Constructs a matrix of threshold candidates (cutpoints) for each feature.
#'   By default, this function attempts to guess if a feature is "binary" 
#'   or "continuous", but you may also specify a \code{feature_types} vector 
#'   with entries in \{"binary","ordinal","continuous"\}.
#'
#' @param X A numeric matrix of size (n x p).
#' @param feature_types A character vector of length p, 
#'   each element one of \{"binary", "ordinal", "continuous"\}.
#'   If \code{NULL}, attempts a simple guess: 
#'   if a column has only \{0,1\}, we call it "binary"; else "continuous".
#' @param num_thresholds Maximum number of thresholds to store per feature (for 
#'   \code{"continuous"} or \code{"ordinal"}).
#' @param bin_continuous Logical. If \code{TRUE}, we pick uniform quantiles 
#'   (excluding min/max). If \code{FALSE}, we attempt to extract unique sorted 
#'   values (excluding extremes).
#'
#' @return A numeric matrix of shape \code{(num_thresholds x p)}, where each 
#'   column is the candidate thresholds for that feature, possibly with 
#'   \code{NA_real_} if fewer than \code{num_thresholds} thresholds exist.
#'
#' @examples
#' X <- matrix(runif(100), nrow = 20) # 20 obs x 5 features
#' feature_types <- c("continuous","binary","continuous","ordinal","continuous")
#' thresh_list <- build_thresholds_list(X, feature_types, num_thresholds=10)
#'
#' @export
build_thresholds_list <- function(
    X,
    feature_types = NULL,
    num_thresholds = 50,
    bin_continuous = TRUE
) {
  if (!is.matrix(X)) {
    stop("X must be a numeric matrix.")
  }
  n <- nrow(X)
  p <- ncol(X)
  
  # If feature_types not given, attempt a naive guess
  if (is.null(feature_types)) {
    feature_types <- vapply(seq_len(p), function(j) {
      vals <- unique(X[, j])
      if (all(vals %in% c(0,1))) {
        "binary"
      } else {
        "continuous"
      }
    }, FUN.VALUE = character(1))
  }
  if (length(feature_types) != p) {
    stop("feature_types must be NULL or length p.")
  }
  
  # Initialize thresholds_list => (num_thresholds x p) of NA
  thresholds_list <- matrix(NA_real_, nrow = num_thresholds, ncol = p)
  
  for (j in seq_len(p)) {
    ftype <- feature_types[j]
    xj <- X[, j]
    
    if (ftype == "binary") {
      # single threshold at 0.5 (or no threshold if you like?)
      thresholds_list[1, j] <- 0.5
      # rest remains NA
    } else if (ftype == "ordinal") {
      # e.g. integer-coded categories => pick midpoints between consecutive 
      # unique sorted values
      vals <- sort(unique(as.numeric(xj)))
      if (length(vals) > 1) {
        cutpoints <- numeric(0)
        for (k in seq_along(vals)[-length(vals)]) {
          midpt <- (vals[k] + vals[k+1]) / 2
          cutpoints <- c(cutpoints, midpt)
        }
        # possibly keep only first num_thresholds
        cutpoints <- head(cutpoints, num_thresholds)
        # fill in thresholds_list
        thresholds_list[seq_along(cutpoints), j] <- cutpoints
      }
    } else if (ftype == "continuous") {
      if (bin_continuous) {
        # pick uniform quantiles
        # e.g. we do num_thresholds + 2 => remove extremes
        if (length(unique(xj)) <= 2) {
          # trivial => do nothing, stays NA
        } else {
          probs <- seq(0, 1, length.out = num_thresholds + 2)
          probs <- probs[-c(1, length(probs))] # remove 0,1 extremes
          cuts  <- as.numeric(quantile(xj, probs = probs, na.rm = TRUE))
          n_cuts <- length(cuts)
          if (n_cuts > num_thresholds) {
            n_cuts <- num_thresholds
          }
          thresholds_list[seq_len(n_cuts), j] <- cuts[seq_len(n_cuts)]
        }
      } else {
        # use unique sorted values, removing min/max
        all_vals <- sort(unique(xj))
        if (length(all_vals) <= 2) {
          # can't do much
        } else {
          all_vals <- all_vals[-c(1, length(all_vals))] # remove extremes
          if (length(all_vals) > num_thresholds) {
            # pick an equally spaced subset
            idx <- round(seq(1, length(all_vals), length.out = num_thresholds))
            all_vals <- all_vals[idx]
          }
          if (length(all_vals) > 0) {
            thresholds_list[seq_len(length(all_vals)), j] <- all_vals
          }
        }
      }
    } else {
      warning(sprintf("Unknown feature type '%s'; leaving column all NA", ftype))
    }
  }
  
  thresholds_list
}
