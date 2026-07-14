### models -- rrBLUP, GBLUP, BayesA, BayesB, Elastic Net, Lasso
# 1. ensure reproducibility
set.seed(2024 + rep_i)

# 2. Parameters
n_folds <- 10  # k-fold CV

# BGLR parameters
bglr_nIter  <- 6000
bglr_burnIn <- 1000
bglr_thin   <- 5

# Paths
data_dir    <- "data/"
output_dir  <- "results/reps/"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 3. Load Packages
library(rrBLUP)
library(BGLR)
library(glmnet)
library(caret)

# 4. Helper Functions

#' Calculate prediction metrics
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

#' rrBLUP model
run_rrBLUP <- function(y_train, geno_train, geno_test) {
  fit  <- mixed.solve(y = y_train, Z = geno_train)
  pred <- as.numeric(geno_test %*% fit$u) + c(fit$beta)
  return(pred)
}

#' GBLUP model
run_GBLUP <- function(y_train, train_names, test_names, G) {
  y_all <- rep(NA, nrow(G))
  names(y_all) <- rownames(G)
  y_all[train_names] <- y_train
  
  dat <- data.frame(gid = names(y_all), y = y_all)
  fit <- kin.blup(data = dat, geno = "gid", pheno = "y", K = G)
  
  return(as.numeric(fit$g[test_names]))
}

#' BayesA model
run_BayesA <- function(y_train, geno_train, geno_test, tmp_prefix) {
  ETA <- list(MRK = list(X = geno_train, model = "BayesA"))
  
  capture.output({
    fit <- BGLR(y = y_train, ETA = ETA,
                nIter = bglr_nIter, burnIn = bglr_burnIn, thin = bglr_thin,
                verbose = FALSE, saveAt = tmp_prefix)
  })
  
  pred <- as.numeric(geno_test %*% fit$ETA$MRK$b) + fit$mu
  
  # Cleanup temp files
  tmp_files <- list.files(dirname(tmp_prefix), 
                          pattern = basename(tmp_prefix), 
                          full.names = TRUE)
  if (length(tmp_files) > 0) file.remove(tmp_files)
  
  return(pred)
}

#' BayesB model
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

#' LASSO model
run_LASSO <- function(y_train, geno_train, geno_test) {
  fit  <- cv.glmnet(x = geno_train, y = y_train, alpha = 1, nfolds = 5)
  pred <- as.numeric(predict(fit, newx = geno_test, s = "lambda.min"))
  return(pred)
}

#' Elastic Net model
run_ElasticNet <- function(y_train, geno_train, geno_test) {
  fit  <- cv.glmnet(x = geno_train, y = y_train, alpha = 0.5, nfolds = 5)
  pred <- as.numeric(predict(fit, newx = geno_test, s = "lambda.min"))
  return(pred)
}

#' Run all models for one fold
run_fold <- function(y, geno, G, train_idx, test_idx, fold_num, pop_name, rep_num) {
  
  y_train    <- y[train_idx]
  y_test     <- y[test_idx]
  geno_train <- geno[train_idx, , drop = FALSE]
  geno_test  <- geno[test_idx, , drop = FALSE]
  train_names <- names(y)[train_idx]
  test_names  <- names(y)[test_idx]
  
  results <- data.frame()
  
  # 1. rrBLUP
  tryCatch({
    pred <- run_rrBLUP(y_train, geno_train, geno_test)
    met  <- calc_metrics(y_test, pred)
    results <- rbind(results, data.frame(
      Population = pop_name, Rep = rep_num, Fold = fold_num, Model = "rrBLUP",
      accuracy = met["accuracy"], rmse = met["rmse"], bias = met["bias"]
    ))
  }, error = function(e) cat("  rrBLUP error:", e$message, "\n"))
  
  # 2. GBLUP
  tryCatch({
    pred <- run_GBLUP(y_train, train_names, test_names, G)
    met  <- calc_metrics(y_test, pred)
    results <- rbind(results, data.frame(
      Population = pop_name, Rep = rep_num, Fold = fold_num, Model = "GBLUP",
      accuracy = met["accuracy"], rmse = met["rmse"], bias = met["bias"]
    ))
  }, error = function(e) cat("  GBLUP error:", e$message, "\n"))
  
  # 3. BayesA
  tryCatch({
    tmp_prefix <- paste0(tempdir(), "/BGLR_BA_", pop_name, "_", rep_num, "_", fold_num, "_")
    pred <- run_BayesA(y_train, geno_train, geno_test, tmp_prefix)
    met  <- calc_metrics(y_test, pred)
    results <- rbind(results, data.frame(
      Population = pop_name, Rep = rep_num, Fold = fold_num, Model = "BayesA",
      accuracy = met["accuracy"], rmse = met["rmse"], bias = met["bias"]
    ))
  }, error = function(e) cat("  BayesA error:", e$message, "\n"))
  
  # 4. BayesB
  tryCatch({
    tmp_prefix <- paste0(tempdir(), "/BGLR_BB_", pop_name, "_", rep_num, "_", fold_num, "_")
    pred <- run_BayesB(y_train, geno_train, geno_test, tmp_prefix)
    met  <- calc_metrics(y_test, pred)
    results <- rbind(results, data.frame(
      Population = pop_name, Rep = rep_num, Fold = fold_num, Model = "BayesB",
      accuracy = met["accuracy"], rmse = met["rmse"], bias = met["bias"]
    ))
  }, error = function(e) cat("  BayesB error:", e$message, "\n"))
  
  # 5. LASSO
  tryCatch({
    pred <- run_LASSO(y_train, geno_train, geno_test)
    met  <- calc_metrics(y_test, pred)
    results <- rbind(results, data.frame(
      Population = pop_name, Rep = rep_num, Fold = fold_num, Model = "LASSO",
      accuracy = met["accuracy"], rmse = met["rmse"], bias = met["bias"]
    ))
  }, error = function(e) cat("  LASSO error:", e$message, "\n"))
  
  # 6. Elastic Net
  tryCatch({
    pred <- run_ElasticNet(y_train, geno_train, geno_test)
    met  <- calc_metrics(y_test, pred)
    results <- rbind(results, data.frame(
      Population = pop_name, Rep = rep_num, Fold = fold_num, Model = "ElasticNet",
      accuracy = met["accuracy"], rmse = met["rmse"], bias = met["bias"]
    ))
  }, error = function(e) cat("  ElasticNet error:", e$message, "\n"))
  
  return(results)
}


# 5. Load and Prepare Data

# Load phenotype files
pheno_B <- read.csv(paste0(data_dir, "pheno_B.csv"), stringsAsFactors = FALSE)
pheno_N <- read.csv(paste0(data_dir, "pheno_N.csv"), stringsAsFactors = FALSE)
pheno_f <- read.csv(paste0(data_dir, "pheno_f.csv"), stringsAsFactors = FALSE)

# Load genotype files
geno_B <- read.csv(paste0(data_dir, "geno_B.csv"), stringsAsFactors = FALSE)
geno_N <- read.csv(paste0(data_dir, "geno_N.csv"), stringsAsFactors = FALSE)
geno_f <- read.csv(paste0(data_dir, "geno_f.csv"), stringsAsFactors = FALSE)

#' Prepare population data
prepare_pop <- function(pheno_df, geno_df) {
  # Genotype matrix
  geno_mat <- as.matrix(geno_df[, -1])
  rownames(geno_mat) <- geno_df[[1]]
  
  # Phenotype vector
  y <- pheno_df[[2]]
  names(y) <- pheno_df[[1]]
  
  # Align
  common <- intersect(names(y), rownames(geno_mat))
  y <- y[common]
  geno_mat <- geno_mat[common, ]
  
  # G matrix
  G <- A.mat(geno_mat)
  
  return(list(y = y, geno = geno_mat, G = G))
}

data_B <- prepare_pop(pheno_B, geno_B)
data_N <- prepare_pop(pheno_N, geno_N)
data_f <- prepare_pop(pheno_f, geno_f)

populations <- list(
  Black = data_B,
  Navy  = data_N,
  Whole = data_f
)

# 6. Run Cross-Validation for This Rep
t_start <- Sys.time()

all_results <- data.frame()

for (pop_name in names(populations)) {
  
  y    <- populations[[pop_name]]$y
  geno <- populations[[pop_name]]$geno
  G    <- populations[[pop_name]]$G
  N    <- length(y)
  
  # Create folds
  folds <- createFolds(1:N, k = n_folds, list = TRUE)
  
  for (k in 1:n_folds) {
    
    cat("  Fold", k, "/", n_folds, "...")
    t_fold <- Sys.time()
    
    test_idx  <- folds[[k]]
    train_idx <- setdiff(1:N, test_idx)
    
    fold_results <- run_fold(y, geno, G, train_idx, test_idx, k, pop_name, rep_i)
    all_results  <- rbind(all_results, fold_results)
    
    cat(" done (", round(difftime(Sys.time(), t_fold, units = "secs"), 1), "s)\n", sep = "")
  }
}

t_end <- Sys.time()


# 7. Save Results
output_file <- paste0(output_dir, "cv_rep_", sprintf("%03d", rep_i), ".csv")
write.csv(all_results, output_file, row.names = FALSE)

