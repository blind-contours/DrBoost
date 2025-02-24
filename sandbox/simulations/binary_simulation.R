#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Disjoint Rule Boosting (logic=TRUE) Simulation for Binary Data
# with THREE new 3-literal rules, each using different logic connectors:
#
#   - 3wayA => (X4=1 & X9=1 & X13=1)        => coefficient = +2.7
#   - 3wayB => ((X2=1 & X7=1) | (X15=1))    => coefficient = +1.8
#   - 3wayC => (X5=1 & (X11=1 | X3=1))      => coefficient = -3.2
#
# Also an intercept = -2.0
# We create a separate validation set for early stopping
# and produce multi-panel black-and-white plots for results.
# ------------------------------------------------------------------

library(mvtnorm)    # for correlated normal => binary
library(ggplot2)
library(patchwork)  # for multi-panel
library(dplyr)
library(tidyr)

# Ensure your "disjoint_rule_boost" is loaded or available:
# source("disjoint_rule_boost.R") or devtools::load_all("...")

set.seed(123)

#========================================================
# 0) Simulation Parameters
#========================================================
p                 <- 20
rho               <- 0.05  # correlation among features
n_vals            <- c(300, 500, 1000, 2000, 4000)
B                 <- 10    # replications
overlap_threshold <- 0.2  # Jaccard overlap threshold for "found" rules

# True logistic coefs for our 3 distinct 3-literal rules:
beta_3wayA <-  2.7
beta_3wayB <-  1.8
beta_3wayC <- -3.2
beta_0     <- -2.0

#========================================================
# 1) Generate correlated binary features
#========================================================
gen_binary_features <- function(n, p, rho=0.3) {
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  Z <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=Sigma)
  X <- (Z > 0)*1
  colnames(X) <- paste0("X", seq_len(p))
  X
}

#========================================================
# 2) True logistic outcome: define coverage for 3wayA,B,C
#========================================================
gen_outcome <- function(X, beta_0, beta_3wayA, beta_3wayB, beta_3wayC) {
  # (A) 3-literal => X4 & X9 & X13
  RA <- as.integer((X[,"X4"]==1) & (X[,"X9"]==1) & (X[,"X13"]==1))
  
  # (B) 3-literal => (X2 & X7) OR (X15)
  # We'll interpret as: (X2=1 & X7=1) OR X15=1
  RB <- integer(nrow(X))
  for(i in seq_len(nrow(X))) {
    RB[i] <- ((X[i,"X2"]==1 && X[i,"X7"]==1) || (X[i,"X15"]==1))
  }
  
  # (C) 3-literal => X5=1 & (X11=1 | X3=1)
  RC <- integer(nrow(X))
  for(i in seq_len(nrow(X))) {
    RC[i] <- ((X[i,"X5"]==1) && ((X[i,"X11"]==1) || (X[i,"X3"]==1)))
  }
  
  eta <- beta_0 +
    beta_3wayA * RA +
    beta_3wayB * RB +
    beta_3wayC * RC
  
  p   <- 1 / (1 + exp(-eta))
  Y   <- rbinom(nrow(X), 1, p)
  
  list(Y=Y, RA=RA, RB=RB, RC=RC, eta=eta)
}

#========================================================
# 3) Overlap + coverage utilities
#========================================================
get_coverage <- function(rule_str, X) {
  dfX <- as.data.frame(X)
  inside <- with(dfX, eval(parse(text=rule_str)))
  as.integer(inside)
}

compute_jaccard <- function(coverA, coverB) {
  intersect_ <- sum(coverA==1 & coverB==1)
  union_     <- sum((coverA==1) | (coverB==1))
  if(union_==0) return(0)
  intersect_/union_
}

#========================================================
# 4) Single replicate => logic=TRUE
#    *Uses separate validation set for early stopping*
#========================================================
simulate_one_run <- function(n, p,
                             beta_0, beta_3wayA, beta_3wayB, beta_3wayC,
                             rho=0.1,
                             max_iter=20,
                             K_twoWay=100,
                             K_threeWay=100,
                             overlap_threshold=0.2) {
  
  # (a) Generate train + validation sets, e.g. 70% / 30%
  n_train <- round(0.7 * n)
  n_val   <- n - n_train
  
  X_full <- gen_binary_features(n, p, rho)
  out    <- gen_outcome(X_full, beta_0, beta_3wayA, beta_3wayB, beta_3wayC)
  Y_full <- out$Y
  
  # split indices
  idx <- sample(seq_len(n), size=n_train, replace=FALSE)
  train_idx <- idx
  val_idx   <- setdiff(seq_len(n), idx)
  
  X_train <- X_full[train_idx, , drop=FALSE]
  Y_train <- Y_full[train_idx]
  X_val   <- X_full[val_idx, , drop=FALSE]
  Y_val   <- Y_full[val_idx]
  
  # (b) Fit Disjoint Rule Boosting
  res <- disjoint_rule_boost(
    X_train      = X_train,
    Y_train      = Y_train,
    X_val        = X_val,
    Y_val        = Y_val,
    outcome_type = "binary",
    logic        = TRUE,
    max_iterations = max_iter,
    K_twoWay     = K_twoWay,
    K_threeWay   = K_threeWay
  )
  final_fit       <- res$final_model
  discovered_rules<- res$rule_list
  num_discovered  <- length(discovered_rules)
  
  if(num_discovered == 0) {
    # No discovered rules
    return(list(
      coef_3wayA=0, coef_3wayB=0, coef_3wayC=0,
      TPR_3wayA=NA, FPR_3wayA=NA,
      TPR_3wayB=NA, FPR_3wayB=NA,
      TPR_3wayC=NA, FPR_3wayC=NA,
      num_rules=0,
      best_overlap_3wayA=0,
      best_overlap_3wayB=0,
      best_overlap_3wayC=0
    ))
  }
  
  # (c) coverage info on training set
  coefs_fit <- coef(final_fit)
  names(coefs_fit) <- gsub("`","", names(coefs_fit))
  
  discovered_info <- lapply(discovered_rules, function(rr) {
    rule_str <- rr$rule_string
    cover_vec<- get_coverage(rule_str, X_train)
    rule_coef<- if(rule_str %in% names(coefs_fit)) coefs_fit[[rule_str]] else NA_real_
    list(cover=cover_vec, coef=rule_coef, rule_str=rule_str)
  })
  
  # (d) Compare discovered coverage vs. true coverage:
  # We'll call them RA,RB,RC on the training set
  RA_train <- out$RA[train_idx]
  RB_train <- out$RB[train_idx]
  RC_train <- out$RC[train_idx]
  
  find_matched_coef <- function(true_cover, dinfo_list, threshold=0.2) {
    best_ov <- 0
    best_cf <- 0
    for(dinf in dinfo_list) {
      ov <- compute_jaccard(true_cover, dinf$cover)
      if(ov > best_ov) {
        best_ov <- ov
        best_cf <- dinf$coef
      }
    }
    if(best_ov < threshold || is.na(best_cf)) {
      return(list(coef=0, best_ov=best_ov))
    } else {
      return(list(coef=best_cf, best_ov=best_ov))
    }
  }
  
  matchA <- find_matched_coef(RA_train, discovered_info, overlap_threshold)
  matchB <- find_matched_coef(RB_train, discovered_info, overlap_threshold)
  matchC <- find_matched_coef(RC_train, discovered_info, overlap_threshold)
  
  best_overlap_3wayA <- matchA$best_ov
  best_overlap_3wayB <- matchB$best_ov
  best_overlap_3wayC <- matchC$best_ov
  
  cA <- matchA$coef
  cB <- matchB$coef
  cC <- matchC$coef
  
  # (e) TPR/FPR => among training
  pred_prob  <- predict(final_fit, type="response")
  pred_class <- as.integer(pred_prob > 0.5)
  
  tpr_fpr_for_rule <- function(true_cover, Y_sub, pred_class_sub) {
    denom_pos <- sum(true_cover==1 & Y_sub==1)
    denom_neg <- sum(true_cover==1 & Y_sub==0)
    tp        <- sum(true_cover==1 & Y_sub==1 & pred_class_sub==1)
    fp        <- sum(true_cover==1 & Y_sub==0 & pred_class_sub==1)
    TPR       <- if(denom_pos>0) tp/denom_pos else NA
    FPR       <- if(denom_neg>0) fp/denom_neg else NA
    c(TPR=TPR, FPR=FPR)
  }
  
  Y_train_sub <- Y_train
  tfA <- tpr_fpr_for_rule(RA_train, Y_train_sub, pred_class)
  tfB <- tpr_fpr_for_rule(RB_train, Y_train_sub, pred_class)
  tfC <- tpr_fpr_for_rule(RC_train, Y_train_sub, pred_class)
  
  list(
    coef_3wayA = cA,
    coef_3wayB = cB,
    coef_3wayC = cC,
    TPR_3wayA  = tfA["TPR"], FPR_3wayA = tfA["FPR"],
    TPR_3wayB  = tfB["TPR"], FPR_3wayB = tfB["FPR"],
    TPR_3wayC  = tfC["TPR"], FPR_3wayC = tfC["FPR"],
    num_rules  = num_discovered,
    best_overlap_3wayA = best_overlap_3wayA,
    best_overlap_3wayB = best_overlap_3wayB,
    best_overlap_3wayC = best_overlap_3wayC
  )
}

#========================================================
# 5) Main loop => store results
#========================================================
sim_results <- list()

for(n_cur in n_vals) {
  cat(sprintf("\n=== Running for n=%d ===\n", n_cur))
  
  # We'll store (A,B,C) coefs in columns
  store_coefs  <- matrix(NA, nrow=B, ncol=3,
                         dimnames=list(NULL, c("coef_3wayA","coef_3wayB","coef_3wayC")))
  # TPR/FPR => 6 columns
  store_tprfpr <- matrix(NA, nrow=B, ncol=6,
                         dimnames=list(NULL, c("TPR_3wayA","FPR_3wayA",
                                               "TPR_3wayB","FPR_3wayB",
                                               "TPR_3wayC","FPR_3wayC")))
  store_numrules <- numeric(B)
  
  # Overlaps => (3 columns: 3wayA, 3wayB, 3wayC)
  store_overlaps <- matrix(0, nrow=B, ncol=3,
                           dimnames=list(NULL,
                                         c("overlap_3wayA","overlap_3wayB","overlap_3wayC")))
  
  for(b in seq_len(B)) {
    outb <- simulate_one_run(
      n           = n_cur,
      p           = p,
      beta_0      = beta_0,
      beta_3wayA  = beta_3wayA,
      beta_3wayB  = beta_3wayB,
      beta_3wayC  = beta_3wayC,
      rho         = rho,
      max_iter    = 10,
      K_twoWay    = 50000,
      K_threeWay  = 50000,
      overlap_threshold = overlap_threshold
    )
    
    # Coeffs
    store_coefs[b,"coef_3wayA"] <- outb$coef_3wayA
    store_coefs[b,"coef_3wayB"] <- outb$coef_3wayB
    store_coefs[b,"coef_3wayC"] <- outb$coef_3wayC
    
    # TPR/FPR
    store_tprfpr[b,"TPR_3wayA"] <- outb$TPR_3wayA
    store_tprfpr[b,"FPR_3wayA"] <- outb$FPR_3wayA
    store_tprfpr[b,"TPR_3wayB"] <- outb$TPR_3wayB
    store_tprfpr[b,"FPR_3wayB"] <- outb$FPR_3wayB
    store_tprfpr[b,"TPR_3wayC"] <- outb$TPR_3wayC
    store_tprfpr[b,"FPR_3wayC"] <- outb$FPR_3wayC
    
    # # rules
    store_numrules[b] <- outb$num_rules
    
    # Overlaps
    store_overlaps[b,"overlap_3wayA"] <- outb$best_overlap_3wayA
    store_overlaps[b,"overlap_3wayB"] <- outb$best_overlap_3wayB
    store_overlaps[b,"overlap_3wayC"] <- outb$best_overlap_3wayC
  }
  
  # Summaries: bias, var, MSE
  true_vals <- c(coef_3wayA=beta_3wayA, coef_3wayB=beta_3wayB, coef_3wayC=beta_3wayC)
  disc_means <- colMeans(store_coefs, na.rm=TRUE)
  bias_vec   <- disc_means - true_vals
  var_vec    <- apply(store_coefs, 2, var, na.rm=TRUE)
  mse_vec    <- colMeans((t(t(store_coefs) - true_vals))^2, na.rm=TRUE)
  
  coef_summary <- data.frame(
    n         = n_cur,
    rule      = c("3wayA","3wayB","3wayC"),
    trueVal   = c(beta_3wayA,beta_3wayB,beta_3wayC),
    meanEst   = disc_means,
    bias      = bias_vec,
    var       = var_vec,
    mse       = mse_vec
  )
  
  # TPR/FPR => average across B
  tprfpr_means <- colMeans(store_tprfpr, na.rm=TRUE)
  tprfpr_df <- data.frame(
    n          = n_cur,
    TPR_3wayA  = tprfpr_means["TPR_3wayA"], FPR_3wayA = tprfpr_means["FPR_3wayA"],
    TPR_3wayB  = tprfpr_means["TPR_3wayB"], FPR_3wayB = tprfpr_means["FPR_3wayB"],
    TPR_3wayC  = tprfpr_means["TPR_3wayC"], FPR_3wayC = tprfpr_means["FPR_3wayC"]
  )
  
  # Overlap => average best overlap
  overlap_avg <- colMeans(store_overlaps, na.rm=TRUE)
  overlap_df  <- data.frame(
    n              = n_cur,
    overlap_3wayA  = overlap_avg["overlap_3wayA"],
    overlap_3wayB  = overlap_avg["overlap_3wayB"],
    overlap_3wayC  = overlap_avg["overlap_3wayC"]
  )
  
  # # discovered rules => mean + sd
  rules_df <- data.frame(
    n              = n_cur,
    mean_num_rules = mean(store_numrules),
    sd_num_rules   = sd(store_numrules)
  )
  
  sim_results[[as.character(n_cur)]] <- list(
    coef_summary = coef_summary,
    tprfpr       = tprfpr_df,
    overlap      = overlap_df,
    numrules     = rules_df
  )
}

#========================================================
# 6) Combine + Print
#========================================================
final_coef_summaries <- do.call(rbind, lapply(sim_results, `[[`, "coef_summary"))
final_tprfpr         <- do.call(rbind, lapply(sim_results, `[[`, "tprfpr"))
final_overlaps       <- do.call(rbind, lapply(sim_results, `[[`, "overlap"))
final_numrules       <- do.call(rbind, lapply(sim_results, `[[`, "numrules"))

cat("\n=== Binary Logic DRB Simulation: 3-literal rules w/ separate val set ===\n")
cat("\nCoefficient Summaries:\n")
print(final_coef_summaries)

cat("\nTPR/FPR Summaries:\n")
print(final_tprfpr)

cat("\nOverlaps:\n")
print(final_overlaps)

cat("\nNumber of Discovered Rules:\n")
print(final_numrules)

#========================================================
# 7) Multi-Panel B/W Figure with consistent linetypes
#========================================================

# We'll define a single factor order => c("3wayA","3wayB","3wayC"),
# then apply a manual scale that assigns the same style to each.

# A) Bias vs. n
df_bias <- final_coef_summaries[, c("n","rule","bias")]
df_bias$n <- as.numeric(as.character(df_bias$n))
df_bias$rule <- factor(df_bias$rule, levels=c("3wayA","3wayB","3wayC"))

p_bias <- ggplot(df_bias, aes(x=n, y=bias, color=rule, linetype=rule, shape=rule)) +
  geom_line() + geom_point() +
  labs(title="Bias of Recovered Coefficients",
       x="Sample Size (n)", y="Bias (Est. - True)") +
  scale_color_manual(values=c("black","black","black")) +
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  scale_shape_manual(values=c(16,17,18)) +
  theme_bw() +
  theme(legend.position="bottom")

# B) TPR vs. n
df_tpr <- final_tprfpr %>%
  select(n, TPR_3wayA, TPR_3wayB, TPR_3wayC)
df_tpr$n <- as.numeric(as.character(df_tpr$n))
df_tpr_long <- pivot_longer(df_tpr,
                            cols=c("TPR_3wayA","TPR_3wayB","TPR_3wayC"),
                            names_to="rule", values_to="TPR"
)
df_tpr_long$rule <- factor(df_tpr_long$rule,
                           levels=c("TPR_3wayA","TPR_3wayB","TPR_3wayC"),
                           labels=c("3wayA","3wayB","3wayC")
)

p_tpr <- ggplot(df_tpr_long, aes(x=n, y=TPR, color=rule, linetype=rule, shape=rule)) +
  geom_line() + geom_point() +
  labs(title="TPR vs. n", x="n", y="True Positive Rate") +
  scale_color_manual(values=c("black","black","black")) +
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  scale_shape_manual(values=c(16,17,18)) +
  theme_bw() +
  theme(legend.position="bottom")

# C) FPR vs. n
df_fpr <- final_tprfpr %>%
  select(n, FPR_3wayA, FPR_3wayB, FPR_3wayC)
df_fpr$n <- as.numeric(as.character(df_fpr$n))
df_fpr_long <- pivot_longer(df_fpr,
                            cols=c("FPR_3wayA","FPR_3wayB","FPR_3wayC"),
                            names_to="rule", values_to="FPR"
)
df_fpr_long$rule <- factor(df_fpr_long$rule,
                           levels=c("FPR_3wayA","FPR_3wayB","FPR_3wayC"),
                           labels=c("3wayA","3wayB","3wayC")
)

p_fpr <- ggplot(df_fpr_long, aes(x=n, y=FPR, color=rule, linetype=rule, shape=rule)) +
  geom_line() + geom_point() +
  labs(title="FPR vs. n", x="n", y="False Positive Rate") +
  scale_color_manual(values=c("black","black","black")) +
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  scale_shape_manual(values=c(16,17,18)) +
  theme_bw() +
  theme(legend.position="bottom")

# D) Overlap vs. n => best overlap
df_ov <- final_overlaps
df_ov$n <- as.numeric(as.character(df_ov$n))
df_ov_long <- pivot_longer(df_ov,
                           cols=c("overlap_3wayA","overlap_3wayB","overlap_3wayC"),
                           names_to="rule", values_to="Overlap"
)
df_ov_long$rule <- factor(df_ov_long$rule,
                          levels=c("overlap_3wayA","overlap_3wayB","overlap_3wayC"),
                          labels=c("3wayA","3wayB","3wayC")
)

p_ov <- ggplot(df_ov_long, aes(x=n, y=Overlap, color=rule, linetype=rule, shape=rule)) +
  geom_line() + geom_point() +
  labs(title="Best Overlap vs. n", x="Sample Size (n)", y="Jaccard Overlap") +
  scale_color_manual(values=c("black","black","black")) +
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  scale_shape_manual(values=c(16,17,18)) +
  theme_bw() +
  theme(legend.position="bottom")

# E) # discovered rules vs. n => single line
df_nr <- final_numrules
df_nr$n <- as.numeric(as.character(df_nr$n))

p_nr <- ggplot(df_nr, aes(x=n, y=mean_num_rules)) +
  geom_line(color="black", linetype="dashed") +
  geom_point(color="black") +
  geom_errorbar(aes(ymin=mean_num_rules - sd_num_rules,
                    ymax=mean_num_rules + sd_num_rules),
                width=0.1, color="black") +
  labs(title="# Discovered Rules vs. n",
       x="n", y="Mean # Rules (+/- 1 SD)") +
  theme_bw()

# Combine subplots into a 2 x 3 layout
p_spacer <- plot_spacer()
combined_plot <- (p_bias + p_tpr + p_fpr) /
  (p_ov + p_nr + p_spacer) +
  plot_annotation(
    title="Binary Logic Simulation",
    subtitle="(3wayA: X4 & X9 & X13, 3wayB: (X2&X7)|(X15), 3wayC: X5 & (X11|X3))"
  )

ggsave("journal_bw_logic_plot.png", combined_plot, width=12, height=8, dpi=300)

cat("\nSaved multi-panel figure: journal_bw_logic_plot.png\n")
cat("\nDone.\n")
