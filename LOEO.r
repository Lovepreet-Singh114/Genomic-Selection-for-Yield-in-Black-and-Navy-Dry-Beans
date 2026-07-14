#!/usr/bin/env Rscript
## ============================================================================ ##
##  Script: LOEO_env.R                                                          ##
##  Leave-One-Environment-Out CV (Environment = Location × Year)               ##
##  Usage: Rscript LOEO_env.R <rep_number>                                     ##
## ============================================================================ ##

# ==============================================================================
# 1. Setup and Arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript LOEO_env.R <rep_number>")
}

rep_i <- as.integer(args[1])
cat("=======================================================\n")
cat("  Leave-One-Environment-Out CV - Repetition", rep_i, "\n")
cat("=======================================================\n\n")

set.seed(2024 + rep_i)

# ==============================================================================
# 2. Parameters
# ==============================================================================

bglr_nIter  <- 6000
bglr_burnIn <- 1000
bglr_thin   <- 5

data_dir    <- "data/"
output_dir  <- "results/env_cv_reps/"

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
  library(tidyverse)
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

blues_env <- read.csv(paste0(data_dir, "blues_environment.csv"), stringsAsFactors = FALSE)
colnames(blues_env) <- c("Genotype", "BLUE", "Level", "Group")

blues_env <- blues_env %>%
  mutate(
    Year = sub(".*_", "", Group),
    Location = sub("_[0-9]+$", "", Group),
    Env = Group  # Environment = Location_Year
  )

# Load genotype data
geno_f <- read.csv(paste0(data_dir, "geno_f.csv"), stringsAsFactors = FALSE)
geno_mat <- as.matrix(geno_f[, -1])
rownames(geno_mat) <- geno_f[[1]]

cat("Genotype matrix:", nrow(geno_mat), "lines x", ncol(geno_mat), "markers\n")

common_genos <- intersect(unique(blues_env$Genotype), rownames(geno_mat))
cat("Common genotypes:", length(common_genos), "\n")

blues_env <- blues_env %>% filter(Genotype %in% common_genos)
geno_mat <- geno_mat[common_genos, ]

all_envs <- sort(unique(blues_env$Env))
cat("Environments:", paste(all_envs, collapse = ", "), "\n")
cat("Number of environments:", length(all_envs), "\n")

# Pre-compute G matrix
cat("\nComputing G matrix...\n")
G <- A.mat(geno_mat)
cat("G matrix computed.\n\n")

rm(geno_f)
gc()

# ==============================================================================
# 6. Run Leave-One-Environment-Out CV
# ==============================================================================

cat("Starting Leave-One-Environment-Out CV for Rep", rep_i, "...\n")
t_start <- Sys.time()

all_results <- data.frame()

for (test_env in all_envs) {
  
  train_envs <- setdiff(all_envs, test_env)
  scenario_name <- paste0("Predict_", test_env)
  
  # Parse location and year from test env
  test_loc  <- sub("_[0-9]+$", "", test_env)
  test_year <- sub(".*_", "", test_env)
  
  cat("\n--- Scenario:", scenario_name, "---\n")
  
  # Training: all other environments (use per-genotype BLUEs averaged across training envs)
  train_data <- blues_env %>% filter(Env != test_env)
  test_data  <- blues_env %>% filter(Env == test_env)
  
  # Average BLUEs across training environments per genotype
  train_avg <- train_data %>%
    group_by(Genotype) %>%
    summarise(y = mean(BLUE, na.rm = TRUE), .groups = "drop")
  
  # Test: BLUEs from the held-out environment
  test_avg <- test_data %>%
    select(Genotype, y = BLUE)
  
  train_genos <- train_avg$Genotype
  test_genos  <- test_avg$Genotype
  
  cat("  Training genotypes:", length(train_genos), "from", length(train_envs), "environments\n")
  cat("  Test genotypes:", length(test_genos), "in", test_env, "\n")
  
  if (length(test_genos) < 10) {
    cat("  Skipping: too few test genotypes\n")
    next
  }
  
  # Prepare named vectors
  y_train <- train_avg$y; names(y_train) <- train_avg$Genotype
  y_test  <- test_avg$y;  names(y_test)  <- test_avg$Genotype
  
  # Prepare genotype matrices
  geno_train <- geno_mat[train_genos, , drop = FALSE]
  geno_test  <- geno_mat[test_genos, , drop = FALSE]
  
  # --- Run all 6 models ---
  
  models <- list(
    list(name = "rrBLUP", fn = function() run_rrBLUP(y_train, geno_train, geno_test)),
    list(name = "GBLUP", fn = function() {
      all_genos_used <- unique(c(train_genos, test_genos))
      G_sub <- G[all_genos_used, all_genos_used]
      run_GBLUP_preG(y_train, train_genos, test_genos, G_sub)
    }),
    list(name = "BayesA", fn = function() {
      tmp_prefix <- paste0(tempdir(), "/BA_LOEO_", test_env, "_")
      run_BayesA(y_train, geno_train, geno_test, tmp_prefix)
    }),
    list(name = "BayesB", fn = function() {
      tmp_prefix <- paste0(tempdir(), "/BB_LOEO_", test_env, "_")
      run_BayesB(y_train, geno_train, geno_test, tmp_prefix)
    }),
    list(name = "LASSO", fn = function() run_LASSO(y_train, geno_train, geno_test)),
    list(name = "ElasticNet", fn = function() run_ElasticNet(y_train, geno_train, geno_test))
  )
  
  for (m in models) {
    tryCatch({
      pred <- m$fn()
      met  <- calc_metrics(y_test, pred)
      all_results <- rbind(all_results, data.frame(
        Scenario = scenario_name, Test_Env = test_env,
        Test_Location = test_loc, Test_Year = test_year,
        Rep = rep_i, Model = m$name,
        n_train = length(train_genos), n_test = length(test_genos),
        accuracy = met["accuracy"], rmse = met["rmse"], bias = met["bias"]))
      cat(" ", m$name, ": r =", round(met["accuracy"], 3), "\n")
    }, error = function(e) cat(" ", m$name, "error:", e$message, "\n"))
  }
  
  gc()
}

# ==============================================================================
# 7. Save Results
# ==============================================================================

t_end <- Sys.time()
cat("\nRep", rep_i, "completed in", round(difftime(t_end, t_start, units = "mins"), 2), "minutes\n")

output_file <- paste0(output_dir, "env_cv_rep_", sprintf("%03d", rep_i), ".csv")
write.csv(all_results, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n")

cat("\n=== Rep", rep_i, "Complete ===\n")