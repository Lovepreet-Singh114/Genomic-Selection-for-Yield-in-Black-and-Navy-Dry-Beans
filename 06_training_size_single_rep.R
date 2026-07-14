#!/usr/bin/env Rscript
## ============================================================================ ##
##                                                                              ##
##  Script: 06_training_size_single_rep.R                                       ##
##                                                                              ##
##  Project: Genomic Selection for Yield Prediction in Dry Beans               ##
##           Comparison of Parametric Models                                    ##
##                                                                              ##
##  Description:                                                                ##
##    Evaluates how training population size affects prediction accuracy.       ##
##    Tests different training set sizes (10%, 20%, ..., 90%) using 10-fold CV. ##
##    Helps determine minimum training population size for acceptable accuracy. ##
##                                                                              ##
##  Usage:                                                                      ##
##    Rscript 06_training_size_single_rep.R <rep_number>                        ##
##                                                                              ##
## ============================================================================ ##

# ==============================================================================
# 1. Setup and Arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 06_training_size_single_rep.R <rep_number>")
}

rep_i <- as.integer(args[1])
cat("=======================================================\n")
cat("  Training Population Size Optimization - Rep", rep_i, "\n")
cat("=======================================================\n\n")

set.seed(2024 + rep_i)

# ==============================================================================
# 2. Parameters
# ==============================================================================

# Training proportions to test
train_proportions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

# BGLR parameters
bglr_nIter  <- 6000
bglr_burnIn <- 1000
bglr_thin   <- 5

# Paths
data_dir    <- "data/"
output_dir  <- "results/training_size_reps/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# 3. Load Packages
# ==============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(rrBLUP)
  library(BGLR)
  library(glmnet)
})
cat("Packages loaded.\n\n")

# ==============================================================================
# 4. Helper Functions
# ==============================================================================

calc_metrics <- function(observed, predicted) {
  valid <- !is.na(observed) & !is.na(predicted)
  obs  <- observed[valid]
  pred <- predicted[valid]
  
  if (length(obs) < 3) {
    return(c(accuracy = NA, rmse = NA, bias = NA))
  }
  
  accuracy <- cor(obs, pred)
  rmse     <- sqrt(mean((obs - pred)^2))
  bias     <- coef(lm(obs ~ pred))[2]
  
  return(c(accuracy = accuracy, rmse = rmse, bias = bias))
}

run_rrBLUP <- function(y_train, geno_train, geno_test) {
  fit  <- mixed.solve(y = y_train, Z = geno_train)
  pred <- as.numeric(geno_test %*% fit$u) + c(fit$beta)
  return(pred)
}

run_GBLUP_preG <- function(y_train, train_names, test_names, G) {
  y_all <- rep(NA, nrow(G))
  names(y_all) <- rownames(G)
  y_all[train_names] <- y_train
  
  dat <- data.frame(gid = names(y_all), y = y_all)
  fit <- kin.blup(data = dat, geno = "gid", pheno = "y", K = G)
  
  return(as.numeric(fit$g[test_names]))
}

run_BayesA <- function(y_train, geno_train, geno_test, tmp_prefix) {
  ETA <- list(MRK = list(X = geno_train, model = "BayesA"))
  
  capture.output({
    fit <- BGLR(y = y_train, ETA = ETA,
                nIter = bglr_nIter, burnIn = bglr_burnIn, thin = bglr_thin,
                verbose = FALSE, saveAt = tmp_prefix)
  })
  
  pred <- as.numeric(geno_test %*% fit$ETA$MRK$b) + fit$mu
  
  tmp_files <- list.files(dirname(tmp_prefix), 
                          pattern = basename(tmp_prefix), 
                          full.names = TRUE)
  if (length(tmp_files) > 0) file.remove(tmp_files)
  
  return(pred)
}

run_BayesB <- function(y_train, geno_train, geno_test, tmp_prefix) {
  ETA <- list(MRK = list(X = geno_train, model = "BayesB"))
  
  capture.output({
    fit <- BGLR(y = y_train, ETA = ETA,
                nIter = bglr_nIter, burnIn = bglr_burnIn, thin = bglr_thin,
                verbose = FALSE, saveAt = tmp_prefix)
  })
  
  pred <- as.numeric(geno_test %*% fit$ETA$MRK$b) + fit$mu
  
  tmp_files <- list.files(dirname(tmp_prefix), 
                          pattern = basename(tmp_prefix), 
                          full.names = TRUE)
  if (length(tmp_files) > 0) file.remove(tmp_files)
  
  return(pred)
}

run_LASSO <- function(y_train, geno_train, geno_test) {
  fit  <- cv.glmnet(x = geno_train, y = y_train, alpha = 1, nfolds = 5)
  pred <- as.numeric(predict(fit, newx = geno_test, s = "lambda.min"))
  return(pred)
}

run_ElasticNet <- function(y_train, geno_train, geno_test) {
  fit  <- cv.glmnet(x = geno_train, y = y_train, alpha = 0.5, nfolds = 5)
  pred <- as.numeric(predict(fit, newx = geno_test, s = "lambda.min"))
  return(pred)
}

# ==============================================================================
# 5. Load and Prepare Data
# ==============================================================================

cat("Loading data...\n")

# Load phenotype and genotype for whole population
pheno_f <- read.csv(paste0(data_dir, "pheno_f.csv"), stringsAsFactors = FALSE)
geno_f  <- read.csv(paste0(data_dir, "geno_f.csv"), stringsAsFactors = FALSE)

# Process genotype matrix
geno_mat <- as.matrix(geno_f[, -1])
rownames(geno_mat) <- geno_f[[1]]

# Process phenotype
y <- pheno_f[[2]]
names(y) <- pheno_f[[1]]

# Align
common <- intersect(names(y), rownames(geno_mat))
y <- y[common]
geno_mat <- geno_mat[common, ]

N <- length(y)
cat("Total genotypes:", N, "\n")
cat("Markers:", ncol(geno_mat), "\n")

# Pre-compute G matrix
cat("Computing G matrix...\n")
G <- A.mat(geno_mat)
cat("G matrix computed.\n\n")

# Clean up
rm(geno_f, pheno_f)
gc()

# ==============================================================================
# 6. Run Training Size Optimization
# ==============================================================================

cat("Starting training size optimization for Rep", rep_i, "...\n")
t_start <- Sys.time()

all_results <- data.frame()

# For each training proportion
for (prop in train_proportions) {
  
  n_train <- round(N * prop)
  n_test  <- N - n_train
  
  cat("\n--- Training proportion:", prop * 100, "% (n=", n_train, ") ---\n", sep = "")
  
  # Random split
  train_idx <- sample(1:N, n_train)
  test_idx  <- setdiff(1:N, train_idx)
  
  y_train     <- y[train_idx]
  y_test      <- y[test_idx]
  geno_train  <- geno_mat[train_idx, , drop = FALSE]
  geno_test   <- geno_mat[test_idx, , drop = FALSE]
  train_names <- names(y)[train_idx]
  test_names  <- names(y)[test_idx]
  
  # rrBLUP
  tryCatch({
    pred <- run_rrBLUP(y_train, geno_train, geno_test)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Rep = rep_i,
      Train_Proportion = prop,
      N_Train = n_train,
      N_Test = n_test,
      Model = "rrBLUP",
      accuracy = met["accuracy"],
      rmse = met["rmse"],
      bias = met["bias"]))
    cat("  rrBLUP: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  rrBLUP error:", e$message, "\n"))
  
  # GBLUP
  tryCatch({
    pred <- run_GBLUP_preG(y_train, train_names, test_names, G)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Rep = rep_i,
      Train_Proportion = prop,
      N_Train = n_train,
      N_Test = n_test,
      Model = "GBLUP",
      accuracy = met["accuracy"],
      rmse = met["rmse"],
      bias = met["bias"]))
    cat("  GBLUP: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  GBLUP error:", e$message, "\n"))
  
  # BayesA
  tryCatch({
    tmp_prefix <- paste0(tempdir(), "/BA_", prop, "_")
    pred <- run_BayesA(y_train, geno_train, geno_test, tmp_prefix)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Rep = rep_i,
      Train_Proportion = prop,
      N_Train = n_train,
      N_Test = n_test,
      Model = "BayesA",
      accuracy = met["accuracy"],
      rmse = met["rmse"],
      bias = met["bias"]))
    cat("  BayesA: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  BayesA error:", e$message, "\n"))
  
  # BayesB
  tryCatch({
    tmp_prefix <- paste0(tempdir(), "/BB_", prop, "_")
    pred <- run_BayesB(y_train, geno_train, geno_test, tmp_prefix)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Rep = rep_i,
      Train_Proportion = prop,
      N_Train = n_train,
      N_Test = n_test,
      Model = "BayesB",
      accuracy = met["accuracy"],
      rmse = met["rmse"],
      bias = met["bias"]))
    cat("  BayesB: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  BayesB error:", e$message, "\n"))
  
  # LASSO
  tryCatch({
    pred <- run_LASSO(y_train, geno_train, geno_test)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Rep = rep_i,
      Train_Proportion = prop,
      N_Train = n_train,
      N_Test = n_test,
      Model = "LASSO",
      accuracy = met["accuracy"],
      rmse = met["rmse"],
      bias = met["bias"]))
    cat("  LASSO: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  LASSO error:", e$message, "\n"))
  
  # ElasticNet
  tryCatch({
    pred <- run_ElasticNet(y_train, geno_train, geno_test)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Rep = rep_i,
      Train_Proportion = prop,
      N_Train = n_train,
      N_Test = n_test,
      Model = "ElasticNet",
      accuracy = met["accuracy"],
      rmse = met["rmse"],
      bias = met["bias"]))
    cat("  ElasticNet: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  ElasticNet error:", e$message, "\n"))
  
  gc()
}

# ==============================================================================
# 7. Save Results
# ==============================================================================

t_end <- Sys.time()
cat("\nRep", rep_i, "completed in", round(difftime(t_end, t_start, units = "mins"), 2), "minutes\n")

output_file <- paste0(output_dir, "training_size_rep_", sprintf("%03d", rep_i), ".csv")
write.csv(all_results, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n")

cat("\n=== Rep", rep_i, "Complete ===\n")
