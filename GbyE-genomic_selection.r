## 3 approaches
library(tidyverse)
library(rrBLUP)
library(BGLR)
library(metan)

# =============================================================================
# GET REP ID FROM COMMAND LINE
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
rep_id <- as.integer(args[1])
set.seed(rep_id)

cat("=== Rep:", rep_id, "===\n")

# =============================================================================
# READ DATA
# =============================================================================
blues <- read.csv("blues_environment.csv") %>%
  rename(GEN = Name, YIELD = BLUE, ENV = Group) %>%
  select(GEN, ENV, YIELD) %>%
  mutate(GEN = as.factor(GEN), ENV = as.factor(ENV))

geno_raw <- read.csv("geno_f.csv")
geno_names <- geno_raw[[1]]
M <- as.matrix(geno_raw[, -1])
rownames(M) <- geno_names

M_centered <- M - 1
K <- A.mat(M_centered)

# =============================================================================
# FIT GGE & EXTRACT ENVIRONMENT SCORES
# =============================================================================
gge_fit <- gge(blues, ENV, GEN, YIELD,
               centering = "environment", scaling = FALSE, svp = "environment")

obj <- gge_fit$YIELD

env_scores <- as.data.frame(obj$coordenv[, 1:2])
colnames(env_scores) <- c("E_PC1", "E_PC2")
env_scores$ENV <- obj$labelenv

# Environment kernel
env_mat <- as.matrix(env_scores[, c("E_PC1", "E_PC2")])
rownames(env_mat) <- env_scores$ENV
cos_sim <- env_mat %*% t(env_mat) /
  (sqrt(rowSums(env_mat^2)) %o% sqrt(rowSums(env_mat^2)))

E_kernel <- cos_sim
E_eigen <- eigen(E_kernel)
E_eigen$values[E_eigen$values < 1e-6] <- 1e-6
E_kernel <- E_eigen$vectors %*% diag(E_eigen$values) %*% t(E_eigen$vectors)
rownames(E_kernel) <- env_scores$ENV
colnames(E_kernel) <- env_scores$ENV

# Merge env scores into blues
blues <- blues %>% left_join(env_scores %>% select(ENV, E_PC1, E_PC2), by = "ENV")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
align_geno <- function(blues_sub, K) {
  common <- intersect(unique(blues_sub$GEN), rownames(K))
  blues_sub <- blues_sub %>% filter(GEN %in% common)
  K_sub <- K[common, common]
  list(blues = blues_sub, K = K_sub, geno = common)
}

run_rrblup <- function(y_train, M_train, M_test) {
  fit <- mixed.solve(y = y_train, Z = M_train)
  y_hat <- M_test %*% fit$u + as.numeric(fit$beta)
  as.numeric(y_hat)
}

# =============================================================================
# SUBSAMPLE FRACTION
# =============================================================================
subsample_frac <- 0.8

# =============================================================================
# STRATEGY 0: BASELINE rrBLUP
# =============================================================================
run_baseline_loeo <- function(blues, M_centered, K) {
  envs <- unique(blues$ENV)
  results <- data.frame()
  
  for (test_env in envs) {
    train <- blues %>% filter(ENV != test_env)
    test  <- blues %>% filter(ENV == test_env)
    
    train_means <- train %>%
      group_by(GEN) %>%
      summarise(YIELD = mean(YIELD, na.rm = TRUE), .groups = "drop")
    
    aligned <- align_geno(train_means, K)
    
    # Subsample training genotypes
    n_train <- length(aligned$geno)
    samp_idx <- sort(sample(n_train, round(n_train * subsample_frac)))
    samp_geno <- aligned$geno[samp_idx]
    
    test_geno <- intersect(test$GEN, samp_geno)
    test_sub <- test %>% filter(GEN %in% test_geno)
    
    M_tr <- M_centered[samp_geno, ]
    M_te <- M_centered[test_geno, ]
    y_tr <- aligned$blues$YIELD[match(samp_geno, aligned$blues$GEN)]
    
    pred_rr <- run_rrblup(y_tr, M_tr, M_te)
    obs <- test_sub$YIELD[match(test_geno, test_sub$GEN)]
    acc_rr <- cor(pred_rr, obs, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Strategy = "Baseline",
      Test_Env = test_env,
      Model = "rrBLUP",
      Accuracy = acc_rr
    ))
  }
  results
}

# =============================================================================
# STRATEGY 1: REACTION NORM rrBLUP
# =============================================================================
run_reaction_norm_loeo <- function(blues, M_centered, K) {
  envs <- unique(blues$ENV)
  results <- data.frame()
  
  for (test_env in envs) {
    train <- blues %>% filter(ENV != test_env)
    test  <- blues %>% filter(ENV == test_env)
    
    common <- intersect(unique(train$GEN), rownames(K))
    common <- intersect(common, unique(test$GEN))
    train <- train %>% filter(GEN %in% common)
    test  <- test %>% filter(GEN %in% common)
    
    # Subsample training genotypes
    train_geno <- unique(train$GEN)
    n_train <- length(train_geno)
    samp_geno <- sort(sample(train_geno, round(n_train * subsample_frac)))
    
    train <- train %>% filter(GEN %in% samp_geno)
    
    env_fit <- lm(YIELD ~ E_PC1 + E_PC2, data = train)
    y_resid_train <- residuals(env_fit)
    
    train_resid <- data.frame(GEN = train$GEN, resid = y_resid_train) %>%
      group_by(GEN) %>%
      summarise(resid = mean(resid), .groups = "drop")
    
    geno_order <- train_resid$GEN
    M_tr <- M_centered[geno_order, ]
    
    fit_rr <- mixed.solve(y = train_resid$resid, Z = M_tr)
    
    test_unique <- test %>%
      filter(GEN %in% samp_geno) %>%
      distinct(GEN, .keep_all = TRUE)
    
    env_pred <- predict(env_fit, newdata = test_unique)
    marker_pred <- M_centered[test_unique$GEN, ] %*% fit_rr$u + as.numeric(fit_rr$beta)
    pred_rr <- env_pred + as.numeric(marker_pred)
    obs <- test_unique$YIELD
    acc_rr <- cor(pred_rr, obs, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Strategy = "Reaction Norm",
      Test_Env = test_env,
      Model = "rrBLUP",
      Accuracy = acc_rr
    ))
  }
  results
}

# =============================================================================
# STRATEGY 3: ENVIRONMENT KERNEL (BGLR)
# =============================================================================
run_env_kernel_loeo <- function(blues, K, E_kernel) {
  envs <- unique(blues$ENV)
  results <- data.frame()
  
  for (test_env in envs) {
    # Subsample genotypes
    all_geno <- intersect(unique(blues$GEN), rownames(K))
    n_geno <- length(all_geno)
    samp_geno <- sort(sample(all_geno, round(n_geno * subsample_frac)))
    
    dat <- blues %>%
      filter(GEN %in% samp_geno) %>%
      mutate(y = ifelse(ENV == test_env, NA, YIELD))
    
    n <- nrow(dat)
    
    gen_idx <- match(dat$GEN, rownames(K))
    K_obs <- K[gen_idx, gen_idx]
    
    env_idx <- match(dat$ENV, rownames(E_kernel))
    E_obs <- E_kernel[env_idx, env_idx]
    
    GxE_obs <- K_obs * E_obs
    
    ETA <- list(
      G = list(K = K_obs, model = "RKHS"),
      GxE = list(K = GxE_obs, model = "RKHS")
    )
    
    fit <- BGLR(
      y = dat$y,
      ETA = ETA,
      nIter = 5000,
      burnIn = 1000,
      verbose = FALSE
    )
    
    test_idx <- which(dat$ENV == test_env)
    pred <- fit$yHat[test_idx]
    obs_df <- blues %>% filter(ENV == test_env, GEN %in% samp_geno)
    obs <- obs_df$YIELD[match(dat$GEN[test_idx], obs_df$GEN)]
    
    acc <- cor(pred, obs, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Strategy = "Env Kernel (BGLR)",
      Test_Env = test_env,
      Model = "GBLUP+GxE",
      Accuracy = acc
    ))
    
    unlink(list.files(pattern = "^ETA|^mu|^var"))
  }
  results
}

# =============================================================================
# RUN ALL STRATEGIES
# =============================================================================
cat("\n=== Running Baseline rrBLUP ===\n")
res_baseline <- run_baseline_loeo(blues, M_centered, K)

cat("\n=== Running Reaction Norm rrBLUP ===\n")
res_reaction <- run_reaction_norm_loeo(blues, M_centered, K)

cat("\n=== Running Environment Kernel ===\n")
res_ekernel <- run_env_kernel_loeo(blues, K, E_kernel)

# =============================================================================
# COMBINE AND SAVE PER-REP RESULTS
# =============================================================================
all_results <- bind_rows(res_baseline, res_reaction, res_ekernel)
all_results$Rep <- rep_id

write.csv(all_results, paste0("results/g_e/GGE_GS_rep_", rep_id, ".csv"),
          row.names = FALSE)

cat("\n=== Rep", rep_id, "Complete ===\n")