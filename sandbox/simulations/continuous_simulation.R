#!/usr/bin/env Rscript
# ----------------------------------------------------------------
# Disjoint Rule Boosting Simulation (logic=FALSE) in [0,1]^2,
# with 4 known piecewise-constant rectangles.
#
# This version produces black-and-white plots:
#  - Distinct line types for each region
#  - Grayscale color scale
#  - Legend labels reflect the actual bounding-box conditions
#
# Author: (Your Name), Date
# ----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(123)

#===============================================================
# 1) Simulation Setup
#===============================================================
# Regions:
#    R1: 0 <= X1 < 0.4, 0 <= X2 < 0.4
#    R2: 0.4 <= X1 < 0.7, 0 <= X2 < 0.6
#    R3: 0.7 <= X1 < 1.0, 0 <= X2 < 0.3
#    R4: everything else

region_labels <- c(
  "1" = "R1: 0<=X1<0.4,\n     0<=X2<0.4",
  "2" = "R2: 0.4<=X1<0.7,\n     0<=X2<0.6",
  "3" = "R3: 0.7<=X1<1.0,\n     0<=X2<0.3",
  "4" = "R4: all other"
)

assign_region_id <- function(x1, x2) {
  if(x1 < 0.4 && x2 < 0.4) {
    return(1L)
  } else if(x1 >=0.4 && x1<0.7 && x2<0.6) {
    return(2L)
  } else if(x1>=0.7 && x2<0.3) {
    return(3L)
  } else {
    return(4L)
  }
}

# Oracle means:
true_region_means <- c(2.0, -1.5, 3.2, 0.0)

# Simulation parameters:
n_vals            <- c(200, 500, 1000, 2000, 4000)
B                 <- 5
sigma_noise       <- 0.3
max_iter          <- 20
min_obs_pct       <- 0.05
max_splits        <- 20
overlap_threshold <- 0.2

#===============================================================
# 2) Data Gen + Utilities
#===============================================================
gen_continuous_data <- function(n, p=2, sigma=0.5) {
  X <- matrix(runif(n*p), nrow=n, ncol=p)
  colnames(X) <- paste0("X", 1:p)
  
  region_ids <- integer(n)
  for(i in seq_len(n)) {
    x1 <- X[i,1]
    x2 <- X[i,2]
    region_ids[i] <- assign_region_id(x1, x2)
  }
  mu <- true_region_means[region_ids]
  Y  <- mu + rnorm(n, mean=0, sd=sigma)
  list(X=X, Y=Y, region_id=region_ids)
}

compute_jaccard <- function(coverA, coverB) {
  intersect_ <- sum(coverA==1 & coverB==1)
  union_     <- sum( (coverA==1)|(coverB==1) )
  if(union_==0) return(0)
  intersect_/union_
}

get_coverage_box <- function(rule_str, X) {
  dfX <- as.data.frame(X)
  inside <- with(dfX, eval(parse(text=rule_str)))
  as.integer(inside)
}

build_thresholds_list <- function(X, max_splits=100) {
  p <- ncol(X)
  out <- matrix(NA_real_, nrow=max_splits, ncol=p)
  for (j in seq_len(p)) {
    x_j <- X[, j]
    unq <- sort(unique(x_j))
    if (length(unq) <= max_splits) {
      out[seq_along(unq), j] <- unq
    } else {
      idx <- round(seq(1, length(unq), length.out = max_splits))
      out[, j] <- unq[idx]
    }
  }
  out
}

#===============================================================
# 3) Single replicate => DRB => match discovered boxes
#===============================================================
simulate_one_run <- function(n,
                             sigma,
                             max_iter,
                             min_obs_pct,
                             overlap_threshold=0.2) {
  # Generate data
  stuff  <- gen_continuous_data(n=n, p=2, sigma=sigma)
  X      <- stuff$X
  Y      <- stuff$Y
  region <- stuff$region_id
  
  thresholds_list <- build_thresholds_list(X, max_splits = max_splits)
  
  # Fit DRB (logic=FALSE => bounding boxes)
  # Note: we do not separate validation set => overfit risk
  res <- disjoint_rule_boost(
    X_train         = X,
    Y_train         = Y,
    X_val           = X,
    Y_val           = Y,
    outcome_type    = "continuous",
    logic           = FALSE,
    thresholds_list = thresholds_list,
    featureNames    = colnames(X),
    max_iterations  = max_iter,
    min_obs_pct     = min_obs_pct,
    patience        = 3, 
    K1 = 1000, 
    K2 = 2000, 
    K3 = 3000
  )
  final_fit       <- res$final_model
  discovered_rules<- res$rule_list
  num_discovered  <- length(discovered_rules)
  if(num_discovered == 0) {
    return(list(
      matched_coefs = rep(0,4),
      coverage_best = rep(0,4),
      num_rules     = 0
    ))
  }
  
  coefs_fit <- coef(final_fit)
  names(coefs_fit) <- gsub("`", "", names(coefs_fit))
  
  discovered_info <- lapply(discovered_rules, function(rr) {
    rule_str  <- rr$rule_string
    coverageV <- get_coverage_box(rule_str, X)
    region_coef <- if(rule_str %in% names(coefs_fit)) coefs_fit[rule_str] else 0
    list(
      coverage = as.integer(coverageV),
      coef     = as.numeric(region_coef),
      rule_str = rule_str
    )
  })
  
  # Match each true region j=1..4 => best Jaccard overlap
  coverage_j_list <- lapply(1:4, function(j) as.integer(region==j))
  find_best_box_for_region <- function(true_cover, disc_info, overlap_thresh) {
    best_ov <- 0
    best_cf <- 0
    for(di in disc_info) {
      ov <- compute_jaccard(true_cover, di$coverage)
      if(ov>best_ov) {
        best_ov <- ov
        best_cf <- di$coef
      }
    }
    if(best_ov < overlap_thresh) {
      return(c(best_cf=0, best_ov=best_ov))
    } else {
      return(c(best_cf=best_cf, best_ov=best_ov))
    }
  }
  
  matched_coefs  <- numeric(4)
  coverage_bests <- numeric(4)
  for(j in 1:4) {
    jj <- find_best_box_for_region(coverage_j_list[[j]], discovered_info, overlap_threshold)
    matched_coefs[j]  <- jj["best_cf"]
    coverage_bests[j] <- jj["best_ov"]
  }
  
  list(
    matched_coefs = matched_coefs,
    coverage_best = coverage_bests,
    num_rules     = num_discovered
  )
}

#===============================================================
# 4) Main loop => multiple n, multiple replicates
#===============================================================
n_vals <- c(200, 500, 1000, 2000)
B      <- 5

sim_results <- list()

for(n_cur in n_vals) {
  cat(sprintf("\n=== Running for n=%d ===\n", n_cur))
  store_coefs    <- matrix(0, nrow=B, ncol=4)
  store_overlaps <- matrix(0, nrow=B, ncol=4)
  store_numrules <- numeric(B)
  
  for(b in seq_len(B)) {
    outb <- simulate_one_run(
      n           = n_cur,
      sigma       = sigma_noise,
      max_iter    = max_iter,
      min_obs_pct = min_obs_pct,
      overlap_threshold = overlap_threshold
    )
    store_coefs[b,]    <- outb$matched_coefs
    store_overlaps[b,] <- outb$coverage_best
    store_numrules[b]  <- outb$num_rules
  }
  
  # Bias–Variance–MSE decomposition
  true_means <- c(2.0, -1.5, 3.2, 0.0)
  disc_mean  <- colMeans(store_coefs)
  bias_vec   <- disc_mean - true_means
  var_vec    <- apply(store_coefs, 2, var)
  mse_vec    <- var_vec + bias_vec^2
  
  region_summary <- data.frame(
    n         = n_cur,
    region_id = 1:4,
    trueMean  = true_means,
    meanEst   = disc_mean,
    bias      = bias_vec,
    var       = var_vec,
    mse       = mse_vec
  )
  
  overlap_means <- colMeans(store_overlaps)
  nr_mean <- mean(store_numrules)
  nr_sd   <- sd(store_numrules)
  
  sim_results[[as.character(n_cur)]] <- list(
    region_summary = region_summary,
    overlap_means  = data.frame(
      n           = n_cur,
      region_id   = 1:4,
      mean_overlap= overlap_means
    ),
    rules_summary  = data.frame(
      n              = n_cur,
      mean_num_rules = nr_mean,
      sd_num_rules   = nr_sd
    )
  )
}

#===============================================================
# 5) Combine & Print
#===============================================================
final_region_summaries <- do.call(rbind, lapply(sim_results, `[[`, "region_summary"))
final_overlaps         <- do.call(rbind, lapply(sim_results, `[[`, "overlap_means"))
final_rules            <- do.call(rbind, lapply(sim_results, `[[`, "rules_summary"))

cat("\n--- Continuous Partition Simulation: Summaries ---\n")
cat("\nRegion-level coefficient estimates:\n")
print(final_region_summaries)

cat("\nMean overlap coverage:\n")
print(final_overlaps)

cat("\nNumber of discovered rules:\n")
print(final_rules)

#===============================================================
# 6) Multi-panel figure in black & white w/ dashed lines
#===============================================================
library(tidyr)

# region_labels => define your bounding boxes for the legend
region_labels <- c(
  "1" = "R1: 0<=X1<0.4,\n     0<=X2<0.4",
  "2" = "R2: 0.4<=X1<0.7,\n     0<=X2<0.6",
  "3" = "R3: 0.7<=X1<1.0,\n     0<=X2<0.3",
  "4" = "R4: all other"
)

#######################
# (a) Bias vs. n
#######################
df_bias <- final_region_summaries[, c("n","region_id","bias")]
df_bias$n <- as.numeric(as.character(df_bias$n))

p_bias <- ggplot(df_bias, aes(x=n, y=bias,
                              group=factor(region_id),
                              color=factor(region_id),
                              linetype=factor(region_id))) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray50") +  # Add zero line
  geom_line() + geom_point() +
  labs(title="Bias of Matched Coefficients by Region",
       x="Sample Size (n)",
       y="Bias (Est. - True)") +
  scale_color_grey(
    start=0.0, end=0.6,
    name="True Region",
    breaks=c("1","2","3","4"),
    labels=region_labels
  ) +
  scale_linetype_discrete(
    name="True Region",
    breaks=c("1","2","3","4"),
    labels=region_labels
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

#######################
# (b) Coverage vs. n
#######################
df_ov <- final_overlaps
df_ov$n <- as.numeric(as.character(df_ov$n))

p_cov <- ggplot(df_ov, aes(x=n, y=mean_overlap,
                           group=factor(region_id),
                           color=factor(region_id),
                           linetype=factor(region_id))) +
  geom_line() + geom_point() +
  labs(title="Coverage (Jaccard) by Region",
       x="Sample Size (n)",
       y="Mean Overlap") +
  scale_color_grey(
    start=0.0, end=0.6,
    name="True Region",
    breaks=c("1","2","3","4"),
    labels=region_labels
  ) +
  scale_linetype_discrete(
    name="True Region",
    breaks=c("1","2","3","4"),
    labels=region_labels
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

#######################
# (c) # discovered rules vs. n
#######################
df_nr <- final_rules
df_nr$n <- as.numeric(as.character(df_nr$n))

# Single line => we just do black color
p_nr <- ggplot(df_nr, aes(x=n, y=mean_num_rules)) +
  geom_line(color="black", linetype="dashed") +
  geom_point(color="black") +
  geom_errorbar(aes(ymin=mean_num_rules - sd_num_rules,
                    ymax=mean_num_rules + sd_num_rules),
                width=0.1, color="black") +
  labs(title="# of Discovered Boxes vs. n",
       x="Sample Size (n)",
       y="Mean # Boxes (+/- 1 SD)") +
  theme_bw()

#######################
# (d) MSE vs. n => average across regions
#######################
df_mse <- final_region_summaries[, c("n","region_id","mse")]
df_mse$n <- as.numeric(as.character(df_mse$n))
df_mse_agg <- df_mse %>%
  group_by(n) %>%
  summarize(mean_mse = mean(mse))

p_mse <- ggplot(df_mse_agg, aes(x=n, y=mean_mse)) +
  geom_line(color="black", linetype="dashed") +
  geom_point(color="black") +
  labs(title="Mean MSE (across all 4 Regions)",
       x="Sample Size (n)",
       y="MSE = Var + Bias^2") +
  theme_bw()

# Combine with patchwork in 2x2 layout
combined_plot <- (p_bias + p_cov) / (p_nr + p_mse) +
  plot_annotation(title="DRB Continuous Simulation",
                  subtitle="Oracle Partition in [0,1]^2, logic=FALSE")

ggsave("journal_bw_plot.png", combined_plot, width=10, height=8, dpi=300)

cat("\nSaved multi-panel figure: journal_bw_plot.png\n")
cat("\nDone.\n")
