################################################################################
# Environmental Ensemble Genomic Selection Models
# Based on Chiaravallotti et al. (2025) - The Plant Genome
# Adapted for MSU Dry Bean Breeding Program
# 
# Three approaches under LOEO cross-validation:
#   1. Singular model  - pool all training environments into one rrBLUP
#   2. Ensemble model  - one rrBLUP per training environment, average predictions
#   3. Optimized ensemble - keep only submodels that improve accuracy
################################################################################

library(rrBLUP)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

# Genotype matrix
geno <- read.csv("geno_f.csv", check.names = FALSE)
rownames(geno) <- geno$ID
M <- as.matrix(geno[, -1])  # 483 x 4986

# BLUEs per environment
blues <- read.csv("blues_environment.csv", stringsAsFactors = FALSE)
colnames(blues) <- c("Name", "BLUE", "Level", "Group")

# Verify environments
envs <- sort(unique(blues$Group))
cat("Environments:", envs, "\n")
cat("Number of environments:", length(envs), "\n")

# Keep only genotypes present in both files
common_ids <- intersect(blues$Name, rownames(M))
cat("Genotypes in BLUEs:", length(unique(blues$Name)), "\n")
cat("Genotypes in marker matrix:", nrow(M), "\n")
cat("Common genotypes:", length(common_ids), "\n")

blues <- blues %>% filter(Name %in% common_ids)
M_common <- M[common_ids, ]

# ── 2. HELPER FUNCTION: rrBLUP PREDICT ───────────────────────────────────────

run_rrblup <- function(train_ids, train_y, pred_ids, marker_matrix) {
  # Fit rrBLUP (mixed.solve) on training data, predict on prediction set
  # train_ids: character vector of training genotype IDs
  # train_y:   named numeric vector of BLUEs for training genotypes
  # pred_ids:  character vector of genotypes to predict
  # marker_matrix: full marker matrix (all genotypes)
  
  M_train <- marker_matrix[train_ids, , drop = FALSE]
  M_pred  <- marker_matrix[pred_ids, , drop = FALSE]
  
  # Solve mixed model
  sol <- mixed.solve(y = train_y, Z = M_train)
  
  # Predict: intercept + M * marker_effects
  gebv <- as.numeric(M_pred %*% sol$u) + sol$beta
  names(gebv) <- pred_ids
  return(gebv)
}

# ── 3. LOEO CROSS-VALIDATION ─────────────────────────────────────────────────

results_list <- list()

for (test_env in envs) {
  
  cat("\n========== Test environment:", test_env, "==========\n")
  
  # Get test set genotypes and their BLUEs
  test_data <- blues %>% filter(Group == test_env)
  test_ids  <- test_data$Name
  test_y    <- setNames(test_data$BLUE, test_data$Name)
  
  # Training environments
  train_envs <- setdiff(envs, test_env)
  
  # ── 3a. SINGULAR MODEL ──────────────────────────────────────────────────
  # Pool all training environments: average BLUEs across environments per genotype
  train_data_all <- blues %>%
    filter(Group %in% train_envs, Name %in% common_ids)
  
  train_blues_pooled <- train_data_all %>%
    group_by(Name) %>%
    summarise(BLUE = mean(BLUE), .groups = "drop")
  
  # Only use genotypes in both training and marker matrix
  train_ids_singular <- intersect(train_blues_pooled$Name, rownames(M_common))
  train_y_singular   <- setNames(
    train_blues_pooled$BLUE[match(train_ids_singular, train_blues_pooled$Name)],
    train_ids_singular
  )
  
  # Predict test genotypes (only those also in marker matrix)
  pred_ids <- intersect(test_ids, rownames(M_common))
  
  if (length(pred_ids) < 5) {
    cat("  Skipping", test_env, "- too few overlapping genotypes\n")
    next
  }
  
  gebv_singular <- run_rrblup(train_ids_singular, train_y_singular, 
                               pred_ids, M_common)
  
  obs_y <- test_y[pred_ids]
  acc_singular <- cor(gebv_singular[pred_ids], obs_y, use = "complete.obs")
  cat("  Singular accuracy:", round(acc_singular, 4), "\n")
  
  # ── 3b. ENSEMBLE MODEL ─────────────────────────────────────────────────
  # Train one rrBLUP per training environment, average predictions
  submodel_preds <- matrix(NA, nrow = length(pred_ids), ncol = length(train_envs),
                           dimnames = list(pred_ids, train_envs))
  
  for (tr_env in train_envs) {
    tr_data <- blues %>% filter(Group == tr_env, Name %in% rownames(M_common))
    tr_ids  <- tr_data$Name
    tr_y    <- setNames(tr_data$BLUE, tr_data$Name)
    
    # Only predict if enough training genotypes
    if (length(tr_ids) < 10) {
      cat("    Skipping submodel", tr_env, "- too few genotypes\n")
      next
    }
    
    gebv_sub <- tryCatch(
      run_rrblup(tr_ids, tr_y, pred_ids, M_common),
      error = function(e) { 
        cat("    Error in submodel", tr_env, ":", e$message, "\n")
        return(NULL) 
      }
    )
    
    if (!is.null(gebv_sub)) {
      submodel_preds[pred_ids, tr_env] <- gebv_sub[pred_ids]
    }
  }
  
  # Average across all submodels (ignore NAs)
  gebv_ensemble <- rowMeans(submodel_preds, na.rm = TRUE)
  acc_ensemble <- cor(gebv_ensemble[pred_ids], obs_y, use = "complete.obs")
  cat("  Ensemble accuracy:", round(acc_ensemble, 4), "\n")
  
  # ── 3c. OPTIMIZED ENSEMBLE ─────────────────────────────────────────────
  # Greedy backward elimination: remove submodels that hurt accuracy
  # Start with full ensemble, drop one at a time, keep if accuracy improves
  
  valid_submodels <- colnames(submodel_preds)[
    colSums(!is.na(submodel_preds)) > 0
  ]
  
  if (length(valid_submodels) >= 2) {
    
    current_set <- valid_submodels
    current_acc <- acc_ensemble
    improved <- TRUE
    
    while (improved && length(current_set) > 1) {
      improved <- FALSE
      best_acc <- current_acc
      best_drop <- NULL
      
      for (candidate in current_set) {
        trial_set <- setdiff(current_set, candidate)
        trial_pred <- rowMeans(submodel_preds[, trial_set, drop = FALSE], na.rm = TRUE)
        trial_acc <- cor(trial_pred[pred_ids], obs_y, use = "complete.obs")
        
        if (!is.na(trial_acc) && trial_acc > best_acc) {
          best_acc <- trial_acc
          best_drop <- candidate
        }
      }
      
      if (!is.null(best_drop)) {
        current_set <- setdiff(current_set, best_drop)
        current_acc <- best_acc
        improved <- TRUE
        cat("    Dropped", best_drop, "-> accuracy:", round(best_acc, 4), "\n")
      }
    }
    
    gebv_optimized <- rowMeans(submodel_preds[, current_set, drop = FALSE], na.rm = TRUE)
    acc_optimized <- cor(gebv_optimized[pred_ids], obs_y, use = "complete.obs")
    cat("  Optimized ensemble accuracy:", round(acc_optimized, 4), "\n")
    cat("  Retained environments:", paste(current_set, collapse = ", "), "\n")
    
  } else {
    acc_optimized <- acc_ensemble
    current_set <- valid_submodels
    gebv_optimized <- gebv_ensemble
  }
  
  # ── Store results ───────────────────────────────────────────────────────
  results_list[[test_env]] <- data.frame(
    TestEnv = test_env,
    Singular = acc_singular,
    Ensemble = acc_ensemble,
    OptimizedEnsemble = acc_optimized,
    N_test = length(pred_ids),
    N_submodels_full = length(valid_submodels),
    N_submodels_opt = length(current_set),
    Retained_envs = paste(current_set, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

# ── 4. COMPILE RESULTS ───────────────────────────────────────────────────────

results_df <- bind_rows(results_list)
print(results_df %>% select(TestEnv, Singular, Ensemble, OptimizedEnsemble, 
                             N_submodels_opt, Retained_envs))

# Summary stats
cat("\n── Mean accuracy across environments ──\n")
cat("Singular:           ", round(mean(results_df$Singular), 4), "\n")
cat("Ensemble:           ", round(mean(results_df$Ensemble), 4), "\n")
cat("Optimized Ensemble: ", round(mean(results_df$OptimizedEnsemble), 4), "\n")

# ── 5. VISUALIZATION ─────────────────────────────────────────────────────────

# Reshape for plotting
plot_df <- results_df %>%
  select(TestEnv, Singular, Ensemble, OptimizedEnsemble) %>%
  pivot_longer(cols = -TestEnv, names_to = "Model", values_to = "Accuracy") %>%
  mutate(
    Model = factor(Model, levels = c("Singular", "Ensemble", "OptimizedEnsemble"),
                   labels = c("Singular rrBLUP", "Ensemble rrBLUP", 
                              "Optimized Ensemble rrBLUP"))
  )

# Bar plot comparison
p1 <- ggplot(plot_df, aes(x = TestEnv, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("Singular rrBLUP" = "#E74C3C", 
                                "Ensemble rrBLUP" = "#3498DB",
                                "Optimized Ensemble rrBLUP" = "#2ECC71")) +
  labs(title = "Environmental Ensemble GS: LOEO Cross-Validation",
       subtitle = "Singular vs Ensemble vs Optimized Ensemble (rrBLUP)",
       x = "Test Environment", y = "Prediction Accuracy (r)",
       fill = "Approach") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

print(p1)
ggsave("ensemble_loeo_comparison.png", p1, width = 10, height = 6, dpi = 300)

# ── 6. DELTA ACCURACY PLOT ───────────────────────────────────────────────────

delta_df <- results_df %>%
  mutate(
    Delta_Ensemble = Ensemble - Singular,
    Delta_Optimized = OptimizedEnsemble - Singular
  ) %>%
  select(TestEnv, Delta_Ensemble, Delta_Optimized) %>%
  pivot_longer(-TestEnv, names_to = "Comparison", values_to = "Delta") %>%
  mutate(Comparison = factor(Comparison, 
                              levels = c("Delta_Ensemble", "Delta_Optimized"),
                              labels = c("Ensemble - Singular", 
                                         "Optimized - Singular")))

p2 <- ggplot(delta_df, aes(x = TestEnv, y = Delta, fill = Comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("Ensemble - Singular" = "#3498DB",
                                "Optimized - Singular" = "#2ECC71")) +
  labs(title = "Change in Accuracy Relative to Singular Model",
       x = "Test Environment", y = expression(Delta ~ "Accuracy (r)"),
       fill = "") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

print(p2)
ggsave("ensemble_delta_accuracy.png", p2, width = 10, height = 6, dpi = 300)

# ── 7. SAVE RESULTS ──────────────────────────────────────────────────────────

write.csv(results_df, "ensemble_gs_results.csv", row.names = FALSE)
cat("\nResults saved to ensemble_gs_results.csv\n")

# ── 8. INDIVIDUAL SUBMODEL ACCURACY TABLE ────────────────────────────────────
# For deeper analysis: how well does each training env predict each test env

cat("\n── Individual submodel accuracy matrix ──\n")

submodel_acc_matrix <- matrix(NA, nrow = length(envs), ncol = length(envs),
                               dimnames = list(paste0("Train_", envs), 
                                               paste0("Test_", envs)))

for (test_env in envs) {
  test_data <- blues %>% filter(Group == test_env)
  test_ids  <- intersect(test_data$Name, rownames(M_common))
  test_y    <- setNames(test_data$BLUE[match(test_ids, test_data$Name)], test_ids)
  
  for (train_env in setdiff(envs, test_env)) {
    tr_data <- blues %>% filter(Group == train_env, Name %in% rownames(M_common))
    tr_ids  <- tr_data$Name
    tr_y    <- setNames(tr_data$BLUE, tr_data$Name)
    
    pred_ids <- intersect(test_ids, rownames(M_common))
    
    if (length(tr_ids) >= 10 && length(pred_ids) >= 5) {
      gebv <- tryCatch(
        run_rrblup(tr_ids, tr_y, pred_ids, M_common),
        error = function(e) NULL
      )
      if (!is.null(gebv)) {
        submodel_acc_matrix[paste0("Train_", train_env), 
                            paste0("Test_", test_env)] <- 
          cor(gebv[pred_ids], test_y[pred_ids], use = "complete.obs")
      }
    }
  }
}

print(round(submodel_acc_matrix, 3))
write.csv(submodel_acc_matrix, "submodel_accuracy_matrix.csv")
cat("Submodel accuracy matrix saved to submodel_accuracy_matrix.csv\n")