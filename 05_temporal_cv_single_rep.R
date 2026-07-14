#!/usr/bin/env Rscript
## ============================================================================ ##
##                                                                              ##
##  Script: 05_temporal_cv_single_rep.R                                         ##
##                                                                              ##
##  Project: Genomic Selection for Yield Prediction in Dry Beans               ##
##           Comparison of Parametric Models                                    ##
##                                                                              ##
##  Description:                                                                ##
##    Runs temporal cross-validation (Leave-One-Year-Out) for all 6            ##
##    parametric models. Evaluates forward prediction scenarios:               ##
##      - Train 2021,2022,2023 → Predict 2024                                  ##
##      - Train 2021,2022,2024 → Predict 2023                                  ##
##      - Train 2021,2023,2024 → Predict 2022                                  ##
##      - Train 2022,2023,2024 → Predict 2021                                  ##
##                                                                              ##
##  Note: This analysis uses genotypes present in both training and test       ##
##  years to enable prediction. Each rep uses different random sampling        ##
##  when genotypes have multiple year observations.                            ##
##                                                                              ##
##  Usage:                                                                      ##
##    Rscript 05_temporal_cv_single_rep.R <rep_number>                          ##
##                                                                              ##
## ============================================================================ ##

# ==============================================================================
# 1. Setup and Arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 05_temporal_cv_single_rep.R <rep_number>")
}

rep_i <- as.integer(args[1])
cat("=======================================================\n")
cat("  Temporal CV (Leave-One-Year-Out) - Repetition", rep_i, "\n")
cat("=======================================================\n\n")

set.seed(2024 + rep_i)

# ==============================================================================
# 2. Parameters
# ==============================================================================

# BGLR parameters
bglr_nIter  <- 6000
bglr_burnIn <- 1000
bglr_thin   <- 5

# Years
all_years <- c("2021", "2022", "2023", "2024")

# Paths
data_dir    <- "data/"
output_dir  <- "results/temporal_cv_reps/"

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

# Load year-specific BLUEs (wide format)
blues_wide <- read.csv(paste0(data_dir, "blues_yearly_wide.csv"), stringsAsFactors = FALSE)

# Load genotype data
geno_f <- read.csv(paste0(data_dir, "geno_f.csv"), stringsAsFactors = FALSE)

# Process genotype matrix
geno_mat <- as.matrix(geno_f[, -1])
rownames(geno_mat) <- geno_f[[1]]

cat("Genotype matrix:", nrow(geno_mat), "lines x", ncol(geno_mat), "markers\n")

# Process BLUEs - get genotype names
blues_names <- blues_wide[[1]]

# Check column names for year BLUEs
cat("BLUE columns:", names(blues_wide), "\n")

# Find common genotypes between BLUEs and genotypes
common_genos <- intersect(blues_names, rownames(geno_mat))
cat("Common genotypes:", length(common_genos), "\n")

# Subset to common genotypes
blues_wide <- blues_wide[blues_wide[[1]] %in% common_genos, ]
geno_mat <- geno_mat[common_genos, ]

# Create named vectors for each year's BLUEs
# Assuming columns are: Name, BLUE_2021, BLUE_2022, BLUE_2023, BLUE_2024
y_2021 <- blues_wide[[2]]; names(y_2021) <- blues_wide[[1]]
y_2022 <- blues_wide[[3]]; names(y_2022) <- blues_wide[[1]]
y_2023 <- blues_wide[[4]]; names(y_2023) <- blues_wide[[1]]
y_2024 <- blues_wide[[5]]; names(y_2024) <- blues_wide[[1]]

# Store in list
y_by_year <- list(
  "2021" = y_2021,
  "2022" = y_2022,
  "2023" = y_2023,
  "2024" = y_2024
)

# Count genotypes per year (non-NA)
for (yr in all_years) {
  n_geno <- sum(!is.na(y_by_year[[yr]]))
  cat("  Year", yr, ":", n_geno, "genotypes with data\n")
}

# Pre-compute G matrix for all genotypes
cat("\nComputing G matrix...\n")
G <- A.mat(geno_mat)
cat("G matrix computed.\n\n")

# Clean up
rm(geno_f, blues_wide)
gc()

# ==============================================================================
# 6. Run Temporal Cross-Validation
# ==============================================================================

cat("Starting temporal CV for Rep", rep_i, "...\n")
t_start <- Sys.time()

all_results <- data.frame()

# Define leave-one-year-out scenarios
scenarios <- list(
  list(train_years = c("2021", "2022", "2023"), test_year = "2024"),
  list(train_years = c("2021", "2022", "2024"), test_year = "2023"),
  list(train_years = c("2021", "2023", "2024"), test_year = "2022"),
  list(train_years = c("2022", "2023", "2024"), test_year = "2021")
)

for (scen in scenarios) {
  
  train_years <- scen$train_years
  test_year   <- scen$test_year
  
  scenario_name <- paste0("Train_", paste(train_years, collapse = "_"), "_Predict_", test_year)
  cat("\n--- Scenario:", scenario_name, "---\n")
  
  # Get test set: genotypes with data in test year
  y_test_full <- y_by_year[[test_year]]
  test_genos <- names(y_test_full)[!is.na(y_test_full)]
  
  cat("  Test genotypes (", test_year, "):", length(test_genos), "\n")
  
  if (length(test_genos) < 10) {
    cat("  Skipping: too few test genotypes\n")
    next
  }
  
  # Get training set: genotypes with data in ANY training year
  # For each genotype, average BLUEs across training years where data exists
  train_genos <- c()
  y_train_list <- list()
  
  for (geno in common_genos) {
    # Get BLUEs for this genotype across training years
    blues_train <- sapply(train_years, function(yr) y_by_year[[yr]][geno])
    
    # If at least one training year has data
    if (sum(!is.na(blues_train)) > 0) {
      train_genos <- c(train_genos, geno)
      # Use mean across available training years
      y_train_list[[geno]] <- mean(blues_train, na.rm = TRUE)
    }
  }
  
  y_train <- unlist(y_train_list)
  
  cat("  Training genotypes:", length(train_genos), "\n")
  
  # Find genotypes that are in BOTH train and test (for GBLUP)
  # But also keep train-only genotypes for training
  test_in_train <- intersect(test_genos, train_genos)
  cat("  Test genotypes also in training:", length(test_in_train), "\n")
  
  # For prediction, we predict test_genos
  y_test <- y_test_full[test_genos]
  
  # Prepare genotype matrices
  geno_train <- geno_mat[train_genos, , drop = FALSE]
  geno_test  <- geno_mat[test_genos, , drop = FALSE]
  
  # rrBLUP
  tryCatch({
    pred <- run_rrBLUP(y_train, geno_train, geno_test)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Scenario = scenario_name, 
      Test_Year = test_year,
      Rep = rep_i, 
      Model = "rrBLUP",
      n_train = length(train_genos),
      n_test = length(test_genos),
      accuracy = met["accuracy"], 
      rmse = met["rmse"], 
      bias = met["bias"]))
    cat("  rrBLUP: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  rrBLUP error:", e$message, "\n"))
  
  # GBLUP
  tryCatch({
    # Use subset of G matrix for genotypes in train or test
    all_genos_used <- unique(c(train_genos, test_genos))
    G_sub <- G[all_genos_used, all_genos_used]
    
    pred <- run_GBLUP_preG(y_train, train_genos, test_genos, G_sub)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Scenario = scenario_name, 
      Test_Year = test_year,
      Rep = rep_i, 
      Model = "GBLUP",
      n_train = length(train_genos),
      n_test = length(test_genos),
      accuracy = met["accuracy"], 
      rmse = met["rmse"], 
      bias = met["bias"]))
    cat("  GBLUP: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  GBLUP error:", e$message, "\n"))
  
  # BayesA
  tryCatch({
    tmp_prefix <- paste0(tempdir(), "/BA_", test_year, "_")
    pred <- run_BayesA(y_train, geno_train, geno_test, tmp_prefix)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Scenario = scenario_name, 
      Test_Year = test_year,
      Rep = rep_i, 
      Model = "BayesA",
      n_train = length(train_genos),
      n_test = length(test_genos),
      accuracy = met["accuracy"], 
      rmse = met["rmse"], 
      bias = met["bias"]))
    cat("  BayesA: r =", round(met["accuracy"], 3), "\n")
  }, error = function(e) cat("  BayesA error:", e$message, "\n"))
  
  # BayesB
  tryCatch({
    tmp_prefix <- paste0(tempdir(), "/BB_", test_year, "_")
    pred <- run_BayesB(y_train, geno_train, geno_test, tmp_prefix)
    met  <- calc_metrics(y_test, pred)
    all_results <- rbind(all_results, data.frame(
      Scenario = scenario_name, 
      Test_Year = test_year,
      Rep = rep_i, 
      Model = "BayesB",
      n_train = length(train_genos),
      n_test = length(test_genos),
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
      Scenario = scenario_name, 
      Test_Year = test_year,
      Rep = rep_i, 
      Model = "LASSO",
      n_train = length(train_genos),
      n_test = length(test_genos),
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
      Scenario = scenario_name, 
      Test_Year = test_year,
      Rep = rep_i, 
      Model = "ElasticNet",
      n_train = length(train_genos),
      n_test = length(test_genos),
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

output_file <- paste0(output_dir, "temporal_cv_rep_", sprintf("%03d", rep_i), ".csv")
write.csv(all_results, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n")

cat("\n=== Rep", rep_i, "Complete ===\n")
